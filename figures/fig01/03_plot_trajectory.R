library(data.table)
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(ggnewscale)
library(ggpubr)
library(viridis)
library(mgcv)
library(tidyr)
library(patchwork)
source(file.path("analysis", "_shared", "cluster_labels.R"))

# ── Helper functions for GAM smoothing ────────────────────────────────────────
make_bin_fill <- function(df_fill, flag_col, n_per_bin = 50) {
  
  score_types <- df_fill %>% distinct(score_type) %>% arrange(score_type) %>% pull(score_type)
  fits <- df_fill %>%
    arrange(score_type) %>%
    group_by(score_type) %>%
    group_map(~ mgcv::gam(score ~ s(mid, bs = "cs"), data = .x, method = "REML")) %>%
    setNames(score_types)

  filtered <- df_fill %>%
    filter(.data[[flag_col]]) %>%
    filter(!is.na(xmin), !is.na(xmax), score_type != "")
  
  if (nrow(filtered) == 0) {
    return(tibble(score_type=character(), bin_id=integer(), mid=numeric(),
                  score_smooth=numeric(), ymin=numeric(), ymax=numeric()))
  }
  
  filtered %>%
    arrange(score_type) %>%
    group_by(score_type, bin_id) %>%
    group_modify(~ {
      if (nrow(.x) == 0) return(tibble())
      xs <- seq(.x$xmin[1], .x$xmax[1], length.out = n_per_bin)
      ys <- as.numeric(predict(fits[[.y$score_type]], newdata = data.frame(mid = xs)))
      tibble(mid = xs, score_smooth = ys, ymin = pmin(0, ys), ymax = pmax(0, ys))
    }) %>%
    ungroup()
}

make_line <- function(df_fill, n = 400) {
  
  score_types <- df_fill %>% distinct(score_type) %>% arrange(score_type) %>% pull(score_type)
  
  fits <- df_fill %>%
    filter(!is.na(mid)) %>%
    arrange(score_type) %>%
    group_by(score_type) %>%
    group_map(~ mgcv::gam(score ~ s(mid, bs = "cs"), data = .x, method = "REML")) %>%
    setNames(score_types)
  
  df_fill %>%
    filter(!is.na(xmin), !is.na(xmax)) %>%
    group_by(score_type) %>%
    summarise(xmin = min(xmin), xmax = max(xmax), .groups = "drop") %>%
    arrange(score_type) %>%
    group_by(score_type) %>%
    group_modify(~ {
      xs <- seq(.x$xmin[1], .x$xmax[1], length.out = n)
      ys <- as.numeric(predict(fits[[.y$score_type]], newdata = data.frame(mid = xs)))
      tibble(mid = xs, score_smooth = ys)
    }) %>%
    ungroup()
}

make_bin_separators <- function(df_rib) {
  df_rib %>%
    group_by(score_type, bin_id) %>%
    summarise(x = max(mid), ymin = min(ymin), ymax = max(ymax), .groups = "drop")
}


# ── Single cell type plot ──────────────────────────────────────────────────────
#
# @param scDRS_long      long-format scDRS scores (from .RData)
# @param bin50_scDRS     bin-level enrichment FDR (from .RData)
# @param tt_res          paired t-test results per bin (from .RData)
# @param range_fdr_comp_global  global range for colour scale (from .RData)
# @param edges           data.frame with bin_id / xmin / xmax
# @param metadata        cell metadata with pseudotime and cluster columns
# @param score1_label    score_type string for score1 (e.g. "GWAS_sex")
# @param score2_label    score_type string for score2 (e.g. "RVAS_sex_fdr001")
# @param score1_name     display name for score1 (e.g. "GWAS")
# @param score2_name     display name for score2 (e.g. "RVAS")
# @param mclust_col      column name for cluster ID (e.g. "mclust25")
# @param mclust_fix      named vector for cluster ID remapping (e.g. c("20"="22"))
# @param my_labels       named vector mapping cluster ID to cell type label
# @param legend_position "none", "bottom", "left", etc.
#
# @return a ggplot object
plot_pseudotime_one <- function(scDRS_long, bin50_scDRS, tt_res,
                                range_fdr_comp_global,
                                edges, metadata,
                                score1_label, score2_label,
                                score1_name, score2_name,
                                mclust_col, mclust_fix = NULL,
                                my_labels,
                                legend_position = "none") {
  
  fdr_cutoff     <- 0.05
  nominal_cutoff <- 0.05
  star_cutoff    <- 0.05
  
  # ── Dynamic rib_type labels based on score names ───────────────────────────
  lab_s1_fdr  <- paste0(score1_name, " FDR < 0.05")
  lab_s2_fdr  <- paste0(score2_name, " FDR < 0.05")
  lab_s1_nom  <- paste0(score1_name, " nominal p < 0.05")
  lab_s2_nom  <- paste0(score2_name, " nominal p < 0.05")
  lab_non     <- "Not significant"
  
  rib_levels <- c(lab_non, lab_s1_fdr, lab_s2_fdr, lab_s1_nom, lab_s2_nom)
  rib_colors <- c(
    "Not significant" = "grey85",
    "#2166AC", "#B2182B", "#92C5DE", "#F4A582"
  )
  names(rib_colors) <- rib_levels
  
  # ── df_fill: bin-level mean score + enrichment significance ───────────────
  df_fill <- scDRS_long %>%
    group_by(score_type, bin_id) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    left_join(edges, by = "bin_id") %>%
    left_join(
      bin50_scDRS %>% select(score_type, bin_id, assoc_mcp, assoc_mcp_fdr),
      by = c("score_type", "bin_id")
    ) %>%
    mutate(
      mid          = (xmin + xmax) / 2,
      sig_fdr      = !is.na(assoc_mcp_fdr) & assoc_mcp_fdr < fdr_cutoff,
      sig_nom      = !is.na(assoc_mcp)     & assoc_mcp     < nominal_cutoff,
      sig_nom_only = sig_nom & !sig_fdr,
      nonsig       = !sig_nom
    )
  message("sig_fdr by score_type:")
  print(df_fill %>% group_by(score_type) %>% 
          summarise(n_sig_fdr = sum(sig_fdr), n_sig_nom = sum(sig_nom_only)))
  
  # ── Ribbon and line data ───────────────────────────────────────────────────
  df_rib_fdr <- make_bin_fill(df_fill, "sig_fdr",      n_per_bin = 50) %>%
    mutate(rib_type = case_when(
      score_type == score1_label ~ lab_s1_fdr,
      score_type == score2_label ~ lab_s2_fdr
    ))
  
  df_rib_nom <- make_bin_fill(df_fill, "sig_nom_only", n_per_bin = 50) %>%
    mutate(rib_type = case_when(
      score_type == score1_label ~ lab_s1_nom,
      score_type == score2_label ~ lab_s2_nom
    ))
  
  df_rib_non <- make_bin_fill(df_fill, "nonsig",       n_per_bin = 50) %>%
    mutate(rib_type = lab_non)
  
  df_rib_all <- bind_rows(df_rib_non, df_rib_nom, df_rib_fdr) %>%
    filter(!is.na(rib_type)) %>%          # drop rows where case_when returned NA
    mutate(rib_type = factor(rib_type, levels = rib_levels))
  
  df_line <- make_line(df_fill, n = 400)
  sep_all <- make_bin_separators(df_rib_all)
  
  # ── Strip: paired t-test per bin ──────────────────────────────────────────
  y_max_global <- max(c(df_line$score_smooth, 0), na.rm = TRUE)
  y_min_global <- min(c(df_line$score_smooth, 0), na.rm = TRUE)
  y_range      <- y_max_global - y_min_global
  strip_height <- 0.05 * y_range
  strip_gap    <- 0.05 * y_range
  
  df_strip <- edges %>%
    left_join(tt_res %>% select(bin_id, fdr_bin_compared), by = "bin_id") %>%
    mutate(
      neglog10_fdr_comp = -log10(pmax(fdr_bin_compared, 1e-300)),
      ymin_strip        = y_max_global + strip_gap,
      ymax_strip        = y_max_global + strip_gap + strip_height,
      signif_bin        = fdr_bin_compared < star_cutoff
    )
  
  # ── Lineage arrow ──────────────────────────────────────────────────────────
  top_strip <- max(df_strip$ymax_strip, na.rm = TRUE)
  y_type    <- top_strip + 0.15 * y_range
  
  # Apply cluster ID remapping if provided (e.g. mclust25==20 → 22)
  if (!is.null(mclust_fix)) {
    for (from in names(mclust_fix)) {
      metadata[[mclust_col]][metadata[[mclust_col]] == as.integer(from)] <- as.integer(mclust_fix[[from]])
    }
  }
  
  df_type <- metadata %>%
    dplyr::select(all_of(mclust_col), average_pseudotime) %>%
    mutate(
      cluster = as.character(.data[[mclust_col]]),
      label   = my_labels[cluster]
    ) %>%
    group_by(label) %>%
    summarise(pseudotime = mean(average_pseudotime, na.rm = TRUE), .groups = "drop") %>%
    arrange(pseudotime) %>%
    mutate(label = factor(label, levels = unique(label)))
  
  df_type_arrow <- df_type %>%
    mutate(xend = lead(pseudotime)) %>%
    filter(!is.na(xend))
  
  colors_use <- c(
    "#5A5156", "#6a00a8", "#AA0DFE", "#FE00FA", "#DEA0FD",
    "#B00068", "#C4451C", "#F6222E", "#F8A19F", "#FEAF16",
    "#90AD1C", "#16FF32", "#1CFFCE", "#2ED9FF", "#3283FE", "#325A9B"
  )
  n_labels       <- length(levels(df_type$label))
  my_colors_named <- setNames(
    colors_use[seq_len(min(n_labels, length(colors_use)))],
    levels(df_type$label)
  )
  
  # ── Assemble plot ──────────────────────────────────────────────────────────
  p <- ggplot() +
    # Ribbons
    geom_ribbon(
      data = df_rib_all,
      aes(x = mid, ymin = ymin, ymax = ymax,
          group = interaction(score_type, bin_id),
          fill  = rib_type),
      alpha = 0.55, color = NA
    ) +
    scale_fill_manual(
      name   = "Enrichment (bin-level)",
      values = rib_colors,
      breaks = rib_levels,
      guide  = guide_legend(order = 1)
    ) +
    
    # Bin separators
    geom_segment(
      data = sep_all,
      aes(x = x, xend = x, y = ymin, yend = ymax),
      inherit.aes = FALSE, color = "white", size = 0.3
    ) +
    
    # Smooth lines
    ggnewscale::new_scale_color() +
    geom_line(
      data = df_line,
      aes(x = mid, y = score_smooth, color = score_type),
      size = 0.5
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    scale_color_manual(
      name   = "scDRS score type",
      values = setNames(c("#2166AC", "#B2182B"), c(score1_label, score2_label)),
      labels = setNames(c(score1_name, score2_name), c(score1_label, score2_label)),
      guide  = guide_legend(order = 2)
    ) +
    
    # Strip: score1 vs score2 per-bin t-test
    ggnewscale::new_scale_fill() +
    geom_rect(
      data = df_strip,
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin_strip, ymax = ymax_strip,
          fill = neglog10_fdr_comp),
      inherit.aes = FALSE, color = NA
    ) +
    scale_fill_gradient(
      name   = paste0(score1_name, " vs ", score2_name, "\n-log10(FDR)"),
      low    = "#FFF9C4",
      high   = "#FEAF16",
      limits = range_fdr_comp_global,
      breaks = scales::pretty_breaks(4),
      guide  = guide_colorbar(order = 4)
    ) +
    geom_rect(
      data = df_strip %>% filter(signif_bin),
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin_strip, ymax = ymax_strip),
      inherit.aes = FALSE, fill = NA, color = "black", size = 0.25
    ) +
    
    # Lineage arrow
    ggnewscale::new_scale_color() +
    geom_segment(
      data = df_type_arrow,
      aes(x = pseudotime, y = y_type, xend = xend, yend = y_type),
      color = "black",
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      lineend = "round", size = 0.25
    ) +
    geom_point(
      data = df_type,
      aes(x = pseudotime, y = y_type, color = label),
      size = 1.5
    ) +
    geom_text(
      data = df_type,
      aes(x = pseudotime, y = y_type, label = label, color = label),
      nudge_y = 0.02 * y_range, size = 2, angle = 45, hjust = 0
    ) +
    scale_color_manual(values = my_colors_named, guide = "none") +
    
    # Theme
    scale_y_continuous(
      breaks = seq(-1, 1, by = 0.5),
      expand = expansion(mult = c(0.02, 0.1))
    ) +
    scale_x_continuous(breaks = seq(0, 30, by = 5)) +
    theme_classic() +
    labs(x = "Average pseudotime", y = "Average scDRS score") +
    theme(
      axis.ticks       = element_line(size = 0.5),
      legend.position  = legend_position,
      legend.direction = "vertical",
      legend.box       = "vertical",
      legend.box.just  = "center",
      legend.spacing.x = unit(0.6, "cm"),
      legend.spacing.y = unit(0.3, "cm")
    )
  
  # Apply font sizes
  ggpar(p,
        font.title = c(6, "bold"), font.subtitle = 6, font.caption = 6,
        font.x = 6, font.y = 6, font.xtickslab = 6, font.ytickslab = 6,
        font.legend = 6)
  
}


# ── Main plot function: EN + IN combined ───────────────────────────────────────
#
# @param cfg          one entry from the comparisons config list
# @param resultdir    base path to result directory
# @param en_labels    named vector mapping EN cluster ID to cell type label
# @param in_labels    named vector mapping IN cluster ID to cell type label
# @param en_mclust_fix  named vector for EN cluster remapping (or NULL)
# @param in_mclust_fix  named vector for IN cluster remapping (or NULL)
plot_comparison <- function(cfg, resultdir,
                            en_labels, in_labels,
                            en_mclust_fix = NULL,
                            in_mclust_fix = NULL) {
  
  # Load .RData for this comparison
  load(file.path(
    resultdir, "WangNature/slingshot",
    paste0("scDRS_RVAS_GWAS_norm_score_", cfg$out_label, ".RData")
  ))
  
  # Load bin edges and metadata for EN and IN
  load_edges <- function(cell_type) {
    bin50 <- fread(file.path(
      resultdir, "WangNature/slingshot", cell_type,
      paste0(cell_type, "_pseudotime_edges_50.txt")
    ))
    bv        <- bin50$V1
    bv[1]     <- bv[1]  - 1e-6
    bv[51]    <- bv[51] + 1e-6
    data.frame(
      bin_id = seq_len(length(bv) - 1),
      xmin   = bv[-length(bv)],
      xmax   = bv[-1]
    )
  }
  
  load_metadata <- function(cell_type) {
    fread(file.path(
      resultdir, "WangNature/slingshot", cell_type,
      paste0(cell_type, "_lineage_slingshot_results.csv")
    )) %>% as.data.frame()
  }
  
  edges_EN    <- load_edges("EN")
  edges_IN    <- load_edges("IN")
  metadata_EN <- load_metadata("EN")
  metadata_IN <- load_metadata("IN")
  
  # EN plot (legend hidden, will be taken from IN)
  p_EN <- plot_pseudotime_one(
    scDRS_long          = scDRS_EN_long,
    bin50_scDRS         = bin50_EN_scDRS,
    tt_res              = tt_res_EN,
    range_fdr_comp_global = range_fdr_comp_global,
    edges               = edges_EN,
    metadata            = metadata_EN,
    score1_label        = cfg$score1_label,
    score2_label        = cfg$score2_label,
    score1_name         = cfg$score1_name,
    score2_name         = cfg$score2_name,
    mclust_col          = "mclust25",
    mclust_fix          = en_mclust_fix,
    my_labels           = en_labels,
    legend_position     = "none"
  )
  
  # IN plot (with legend)
  p_IN <- plot_pseudotime_one(
    scDRS_long          = scDRS_IN_long,
    bin50_scDRS         = bin50_IN_scDRS,
    tt_res              = tt_res_IN,
    range_fdr_comp_global = range_fdr_comp_global,
    edges               = edges_IN,
    metadata            = metadata_IN,
    score1_label        = cfg$score1_label,
    score2_label        = cfg$score2_label,
    score1_name         = cfg$score1_name,
    score2_name         = cfg$score2_name,
    mclust_col          = "mclust22",
    mclust_fix          = in_mclust_fix,
    my_labels           = in_labels,
    legend_position     = "left"
  )
  
  # Save EN individually
  ggsave(
    filename = file.path(resultdir, "WangNature/slingshot/EN",
                         paste0("scDRS_RVAS_GWAS_norm_score_avgpseudotime_sexcov_curve_",
                                cfg$out_label, ".png")),
    plot = p_EN, width = 8, height = 8, units = "in"
  )
  
  # Save IN individually
  ggsave(
    filename = file.path(resultdir, "WangNature/slingshot/IN",
                         paste0("scDRS_RVAS_GWAS_norm_score_avgpseudotime_sexcov_curve_",
                                cfg$out_label, ".png")),
    plot = p_IN, width = 8, height = 10, units = "in"
  )
  
  # Save combined EN / IN
  plot_combined <- p_EN / p_IN
  ggsave(
    filename = file.path(resultdir, "WangNature/slingshot",
                         paste0("scDRS_RVAS_GWAS_norm_score_avgpseudotime_sexcov_curve_",
                                cfg$out_label, ".pdf")),
    plot = plot_combined, width = 6, height = 3, units = "in"
  )
  
  message("Saved: ", cfg$out_label)
  invisible(list(EN = p_EN, IN = p_IN))
}


# ── Cell type labels ───────────────────────────────────────────────────────────
en_labels <- get_en_labels()
in_labels <- get_in_labels()

# ── Comparison config list ─────────────────────────────────────────────────────
#
# Each entry defines one pairwise comparison. To add a new comparison,
# simply append a new list() entry here — no function changes needed.
#
# Fields:
#   score1/2_label  : score_type values in scDRS_long and bin50_scDRS
#                     (must match what scale_color_manual uses)
#   score1/2_name   : display names shown in plot legends
#   out_label       : used to locate the .RData file produced by the data script

comparisons <- list(
  
  list(score1_label = "GWAS_sex",       score2_label = "RVAS_sex_fdr001",
       score1_name  = "GWAS",           score2_name  = "RVAS",
       out_label    = "ASC_Pfdr_253genes_GWAS_vs_fdr001"),
  
  list(score1_label = "GWAS_sex",       score2_label = "RVAS_sex_fdr001",
       score1_name  = "GWAS",           score2_name  = "RVAS",
       out_label    = "ASC_Pfdr_253genes_200GWAS_vs_fdr001"),
  
  list(score1_label = "GWAS_sex",       score2_label = "RVAS_sex_fdr01",
       score1_name  = "GWAS",           score2_name  = "RVAS",
       out_label    = "ASC_Pfdr_416genes_GWAS_vs_fdr01"),
  
  list(score1_label = "GWAS_sex",       score2_label = "RVAS_sex_fdr05",
       score1_name  = "GWAS",           score2_name  = "RVAS",
       out_label    = "ASC_Pfdr_696genes_GWAS_vs_fdr05"),
  
  list(score1_label = "GWAS_sex",       score2_label = "RVAS_sex_fdr1",
       score1_name  = "GWAS",           score2_name  = "RVAS",
       out_label    = "ASC_Pfdr_951genes_GWAS_vs_fdr1")
)


# ── Run all comparisons ────────────────────────────────────────────────────────
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"

for (cfg in comparisons) {
  plot_comparison(
    cfg           = cfg,
    resultdir     = resultdir,
    en_labels     = en_labels,
    in_labels     = in_labels,
    en_mclust_fix = c("20" = "22"),   # EN: mclust25==20 → 22
    in_mclust_fix = c("21" = "15")    # IN: mclust22==21 → 15
  )
}
