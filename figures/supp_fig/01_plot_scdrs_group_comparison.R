library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# ── Merge at Group level ──────────────────────────────────────────────────────
compute_group_heatmap_df <- function(df1, df2,
                                     label1 = "group1",
                                     label2 = "group2") {
  stopifnot(label1 != label2)

  merge(
    df1 %>%
      mutate(assoc_mcp_fdr = p.adjust(assoc_mcp, method = "fdr")) %>%
      dplyr::select(group, assoc_mcz, assoc_mcp, assoc_mcp_fdr) %>%
      dplyr::rename(
        z1   = assoc_mcz,
        p1   = assoc_mcp,
        fdr1 = assoc_mcp_fdr
      ),
    df2 %>%
      mutate(assoc_mcp_fdr = p.adjust(assoc_mcp, method = "fdr")) %>%
      dplyr::select(group, assoc_mcz, assoc_mcp, assoc_mcp_fdr) %>%
      dplyr::rename(
        z2   = assoc_mcz,
        p2   = assoc_mcp,
        fdr2 = assoc_mcp_fdr
      ),
    by = "group"
  ) %>%
    mutate(
      label1_name = label1,
      label2_name = label2
    )
}

plot_group_heatmap <- function(df,
                               label1,
                               label2,
                               group_order = NULL,
                               title = NULL,
                               fill_limits = NULL,
                               fill_colors = heatmap_colors) {
  stopifnot(label1 != label2)
  
  if (is.null(group_order)) {
    group_order <- df %>%
      arrange(desc(pmax(z1, z2, na.rm = TRUE))) %>%
      pull(group)
  }
  
  df_long <- df %>%
    mutate(group = factor(group, levels = group_order)) %>%
    pivot_longer(
      cols = c(z1, z2),
      names_to = "score_type",
      values_to = "z_score"
    ) %>%
    mutate(
      score_type = recode(score_type,
                          z1 = label1,
                          z2 = label2),
      fdr_val = ifelse(score_type == label1, fdr1, fdr2),
      sig_label = case_when(
        fdr_val < 0.001 ~ "***",
        fdr_val < 0.01  ~ "**",
        fdr_val < 0.05  ~ "*",
        TRUE            ~ ""
      ),
      score_type = factor(score_type, levels = c(label1, label2))
    )
  
  if (is.null(fill_limits)) {
    zmax <- max(abs(df_long$z_score), na.rm = TRUE)
    fill_limits <- c(-zmax, zmax)
  }
  
  scale_knots <- c(fill_limits[1], -2, -1, 0, 1, 2, 3, fill_limits[2])
  
  ggplot(df_long, aes(x = group, y = score_type, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 3) +
    scale_fill_gradientn(
      colours = fill_colors,
      values = scales::rescale(scale_knots, from = fill_limits),
      limits = fill_limits,
      oob = scales::squish,
      name = "scDRS z-score"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      plot.title = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
}


# ── Build paired EN/IN plot with shared scale ─────────────────────────────────
make_group_pair_plot <- function(df_en_1, df_en_2,
                                 df_in_1, df_in_2,
                                 label1, label2,
                                 title_en, title_in,
                                 group_order = NULL,
                                 fill_colors = heatmap_colors) {
  if (is.null(group_order)) {
    group_order <- c(
      "First_trimester",
      "Second_trimester",
      "Third_trimester",
      "Infancy",
      "Adolescence"
    )
  }
  
  df_group_EN <- compute_group_heatmap_df(df_en_1, df_en_2, label1, label2)
  df_group_IN <- compute_group_heatmap_df(df_in_1, df_in_2, label1, label2)
  
  all_z <- c(
    df_group_EN$z1, df_group_EN$z2,
    df_group_IN$z1, df_group_IN$z2
  )
  common_zmax <- max(abs(all_z), na.rm = TRUE)
  common_limits <- c(-common_zmax, common_zmax)
  
  p_group_EN <- plot_group_heatmap(
    df_group_EN,
    label1 = label1,
    label2 = label2,
    group_order = group_order,
    title = title_en,
    fill_limits = common_limits,
    fill_colors = fill_colors
  )
  
  p_group_IN <- plot_group_heatmap(
    df_group_IN,
    label1 = label1,
    label2 = label2,
    group_order = group_order,
    title = title_in,
    fill_limits = common_limits,
    fill_colors = fill_colors
  )
  
  theme_fig <- theme_classic() +
    theme(legend.position = "top", text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.title = element_text(size = 7), legend.text = element_text(size = 7),
          legend.title = element_text(size = 7), plot.margin = margin(5, 5, 5, 5))
  
  p_group_EN = p_group_EN + theme_fig
  p_group_IN = p_group_IN + theme_fig
  
  p_group_EN / p_group_IN +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "right")
}


# ── custom colors ─────────────────────────────────────────────────────────────
heatmap_colors <- c(
  "#08306B",
  "#6BAED6",
  "#C6DBEF",
  "white",
  "#F9C6A3",
  "#F07020",
  "#E03060",
  "#C0267A"
)

group_order <- c(
  "First_trimester",
  "Second_trimester",
  "Third_trimester",
  "Infancy",
  "Adolescence"
)

# ── Load Group-level scDRS results ────────────────────────────────────────────
trait <- "ASC_Pfdr"
trait_GWAS <- "ASD_Matoba_2020"

group_scDRS_RVAS_EN <- fread(file.path(
  resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/EN",
  paste0(trait, ".scdrs_group.Group")
))
group_scDRS_RVAS_IN <- fread(file.path(
  resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/IN",
  paste0(trait, ".scdrs_group.Group")
))

group_scDRS_GWAS_EN <- fread(file.path(
  resultdir, "WangNature/scDRS/run_GWAS_v3/downstream_analysis/100_genes/EN",
  paste0(trait_GWAS, ".scdrs_group.Group")
))
group_scDRS_GWAS_IN <- fread(file.path(
  resultdir, "WangNature/scDRS/run_GWAS_v3/downstream_analysis/100_genes/IN",
  paste0(trait_GWAS, ".scdrs_group.Group")
))

df = rbind(group_scDRS_RVAS_EN %>% mutate(lineage = "EN", score_type = "RVAS"), 
           group_scDRS_RVAS_IN %>% mutate(lineage = "IN", score_type = "RVAS"),
           group_scDRS_GWAS_EN %>% mutate(lineage = "EN", score_type = "GWAS"), 
           group_scDRS_GWAS_IN%>% mutate(lineage = "IN", score_type = "GWAS"))

fwrite(df, file = file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/",
                            paste0(trait, ".scdrs_group.Group.csv")))


p_rvas_vs_gwas <- make_group_pair_plot(
  df_en_1 = group_scDRS_RVAS_EN,
  df_en_2 = group_scDRS_GWAS_EN,
  df_in_1 = group_scDRS_RVAS_IN,
  df_in_2 = group_scDRS_GWAS_IN,
  label1 = "RVAS",
  label2 = "GWAS",
  title_en = "EN Group-level scDRS: RVAS vs GWAS",
  title_in = "IN Group-level scDRS: RVAS vs GWAS"
)

p_rvas_vs_gwas

ggsave(
  file.path(resultdir, "WangNature/slingshot",
            paste0("scDRS_RVAS_vs_GWAS_Group_heatmap_", trait_GWAS, ".pdf")),
  p_rvas_vs_gwas,
  width = 5,
  height = 2,
  units = "in"
)
