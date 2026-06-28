library(Seurat)
library(Signac)
library(slingshot)
library(tidyverse)
library(scCustomize)
library(dplyr)
library(data.table)
library(scales)
library(circlize)
library(gridGraphics)
library(gridExtra)
library(cowplot)
source(file.path("analysis", "_shared", "cluster_labels.R"))

# ── Paths ──────────────────────────────────────────────────────────────────────
datapath  <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
en_labels <- get_en_labels()
in_labels <- get_in_labels()

# ══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════════════════════

# ── Lineage plot: node color = Z score, asterisk for FDR < 0.05 ───────────────
#
# @param dat_seurat    Seurat object (EN or IN)
# @param sds           SlingshotDataSet object
# @param lineage_scDRS full scDRS result table (both EN and IN rows),
#                      filtered internally by cell_type prefix
# @param z_lim         global z range for symmetric color scale and point size
# @param trait         trait identifier, used to find score column name
# @param gene_number   gene set size, used in column name and file name
# @param version       scDRS run version, used in column name and file name
# @param cell_type     "EN" or "IN"
# @param my_labels     named vector mapping cluster ID to cell type label
# @param resultdir     base path for saving output PDF
#
# @return the FeaturePlot object invisibly (for legend extraction)
plot_lineage <- function(dat_seurat, sds, lineage_scDRS,
                         z_lim,
                         trait, gene_number, version,
                         cell_type, my_labels,
                         resultdir) {
  
  # ── 1. Filter scDRS result for this cell type ──────────────────────────────
  ct_scDRS <- lineage_scDRS %>%
    dplyr::filter(grepl(paste0("^", cell_type), group)) %>%
    mutate(cluster = as.character(cluster))
  
  # ── 2. UMAP feature plot (used only to extract cell colours) ──────────────
  reduction <- paste0("wnn.", cell_type, ".umap2")
  score_col <- paste0("scDRS_", version, "_", gene_number, "_", trait, "_norm_score")
  
  p_scDRS <- FeaturePlot_scCustom(
    dat_seurat,
    features   = score_col,
    reduction  = reduction,
    colors_use = c("#2166AC", "#4393C3", "#92C5DE", "#CADDED",
                   "#FDDBC7",
                   "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    max.cutoff = 5,
    min.cutoff = -3,
    na_cutoff  = NULL,
    raster     = TRUE
  )
  
  p_build <- ggplot_build(p_scDRS)
  df_umap <- p_build$data[[1]]
  
  # ── 3. Draw UMAP scatter (base R) ─────────────────────────────────────────
  plot(df_umap$x, df_umap$y,
       xlab = "", ylab = "",
       col  = scales::alpha(df_umap$colour, 0.3),
       asp  = NA, pch = 16, cex = 0.15,
       axes = FALSE, xaxs = "i", yaxs = "i")
  
  # ── 4. Compute cluster centres from SlingshotDataSet ──────────────────────
  X             <- reducedDim(sds)
  clusterLabels <- slingClusterLabels(sds)
  MST           <- slingMST(sds)
  clusters      <- rownames(MST)
  clusters_chr  <- as.character(clusters)
  
  centers <- t(vapply(clusters, function(clID) {
    w <- clusterLabels[, clID]
    apply(X, 2, weighted.mean, w = w)
  }, numeric(ncol(X))))
  rownames(centers) <- clusters
  
  # ── 5. Map cluster IDs to Z scores and FDR ────────────────────────────────
  z_values   <- ct_scDRS$assoc_mcz[match(clusters_chr, ct_scDRS$cluster)]
  fdr_values <- ct_scDRS$assoc_mcp_fdr[match(clusters_chr, ct_scDRS$cluster)]
  size_val   <- abs(z_values)
  
  # ── 6. Map Z scores to colours (diverging, continuous gradient) ────────────
  # Symmetric color scale centred at 0: blue = negative, red = positive
  z_neg <- z_lim[z_lim < 0]
  z_pos <- z_lim[z_lim > 0]
  
  z_neg_lim <- quantile(z_neg, 0.05, na.rm = TRUE)
  z_pos_lim <- quantile(z_pos, 0.95, na.rm = TRUE)
  
  col_fun <- colorRamp2(
    c(z_lim[1], -1, -0.5, 0, 1.3, 2, 3.5, z_lim[2]),
    c("#08306B",     
      "#6BAED6", 
      "#C6DBEF",   
      "white",
      "#F9C6A3",
      "#F07020",  
      "#E03060", 
      "#C0267A")
  )
  
  pt_cols <- col_fun(ifelse(is.na(z_values), 0, z_values))
  
  # ── 7. Significance: FDR < 0.05 marked with asterisk ─────────────────────
  is_sig    <- !is.na(fdr_values) & fdr_values < 0.05
  pt_border <- "black"
  
  # ── 8. Map |Z| to point sizes ─────────────────────────────────────────────
  cex_range     <- c(0.8, 2.3)
  size_val      <- abs(z_values)
  size_val_clip <- pmin(size_val, max(abs(z_lim)))
  cex_val       <- cex_range[1] +
    (size_val_clip - 0) / max(abs(z_lim)) * diff(cex_range)
  
  # ── 9. Overlay slingshot lineages and cluster points ───────────────────────
  plot(sds,
       type = "lineages", add = TRUE,
       col = "black", cex = 0, lwd = 0.7, asp = 1)
  
  points(centers[, 1], centers[, 2],
         pch = 21, bg = pt_cols, col = pt_border,
         cex = cex_val, lwd = 0.7)
  
  # Three-level significance annotation
  sig_label <- case_when(
    !is.na(fdr_values) & fdr_values < 0.001 ~ "***",
    !is.na(fdr_values) & fdr_values < 0.01  ~ "**",
    !is.na(fdr_values) & fdr_values < 0.05  ~ "*",
    TRUE                                      ~ ""
  )
  
  sig_idx <- sig_label != ""
  if (any(sig_idx)) {
    text(centers[sig_idx, 1], centers[sig_idx, 2],
         labels = sig_label[sig_idx],
         cex    = ifelse(sig_label[sig_idx] == "*", 0.6, 0.5),
         col    = "white",
         font   = 2)
  }
  
  # ── 10. Size legend ────────────────────────────────────────────────────────
  legend_z   <- pretty(z_lim, n = 3)
  legend_z   <- legend_z[legend_z > 0]
  legend_cex <- cex_range[1] + (legend_z - z_lim[1]) / diff(z_lim) * diff(cex_range)
  
  legend_pos <- if (cell_type == "EN") "topleft" else "topright"
  legend(legend_pos,
         legend = paste0("|Z| = ", round(legend_z, 1)),
         pch = 21, pt.cex = legend_cex,
         pt.bg = "grey70", col = "black",
         title = "RVAS enrichment",
         bty = "n", cex = 0.8)
  
  # ── 11. Capture base-R plot as grob and save PDF ──────────────────────────
  gridGraphics::grid.echo()
  slingshot_grob <- grid::grid.grab()
  
  out_dir  <- file.path(resultdir, "WangNature/slingshot", cell_type)
  pdf_size <- if (cell_type == "EN") c(3.8, 3.3) else c(3.5, 3.5)
  
  pdf(file   = file.path(out_dir,
                         paste0("lineage_scDRS_", version, "_norm_score_umap_",
                                gene_number, "genes_pval_", trait, ".pdf")),
      width  = pdf_size[1],
      height = pdf_size[2])
  gridExtra::grid.arrange(slingshot_grob)
  dev.off()
  
  message("Saved lineage plot: ", cell_type, " / ", version, " / ",
          gene_number, " genes / ", trait)
  
  invisible(p_scDRS)
}


# ── Legend: scDRS score colorbar + Z score colorbar + size legend ─────────────
#
# @param p_scDRS     FeaturePlot returned by plot_lineage()
# @param z_lim       global z range (used for symmetric color scale)
# @param version     used in output file name
# @param gene_number used in output file name
# @param trait       used in output file name
# @param resultdir   base path for saving output PDF
# @param cell_type   used in output path
save_lineage_legend <- function(p_scDRS, z_lim,
                                version, gene_number, trait,
                                resultdir, cell_type = "EN") {
  
  z_abs_max <- max(abs(z_lim), na.rm = TRUE)
  
  # ── scDRS normalized score colorbar ─────────────────────────────────────────
  p_scDRS_leg <- p_scDRS +
    labs(color = "scDRS normalized\nscore") +
    guides(color = guide_colorbar(
      title.position = "left", title.hjust = 0.5, title.vjust = 0.5,
      barwidth  = unit(0.35, "cm"),
      barheight = unit(3.2, "cm"),
      ticks     = TRUE
    )) +
    theme_void(base_size = 6) +
    theme(
      legend.position  = "left",
      legend.direction = "vertical",
      legend.title     = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      legend.text      = element_text(size = 6),
      plot.margin      = margin(0, 0, 0, 0)
    )
  legend_scDRS_grob <- cowplot::get_legend(p_scDRS_leg)
  
  # ── RVAS enrichment Z score colorbar (diverging, continuous) ─────────────────
  # In save_lineage_legend(): replace the RVAS enrichment colorbar section
    
  z_breaks_leg <- c(z_lim[1], -1, -0.5, 0, 1.3, 2, 3.5, z_lim[2])
  z_colors_leg <- c("#08306B","#6BAED6", "#C6DBEF", "white",
                    "#F9C6A3","#F07020",  "#E03060", "#C0267A")
  
  p_zscore_leg <- ggplot(data.frame(y = seq(z_lim[1], z_lim[2], length.out = 100)),
                         aes(x = 1, y = y, fill = y)) +
    geom_raster() +
    scale_fill_gradientn(
      colours = z_colors_leg,
      values  = scales::rescale(z_breaks_leg, to = c(0, 1)),
      limits  = c(z_lim[1], z_lim[2]),
      breaks = c(-2,0,2,4),
      name    = "RVAS enrichment\nZ score"
    ) +
    guides(fill = guide_colorbar(
      title.position = "left", title.hjust = 0.5, title.vjust = 0.5,
      barwidth  = unit(0.35, "cm"),
      barheight = unit(3.2, "cm"),
      ticks     = TRUE,
      draw.ulim = TRUE,
      draw.llim = TRUE
    )) +
    theme_void(base_size = 6) +
    theme(
      legend.position  = "left",
      legend.direction = "vertical",
      legend.title     = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 6),
      legend.text      = element_text(size = 5),
      plot.margin      = margin(0, 0, 0, 0)
    )
  legend_zscore_grob <- cowplot::get_legend(p_zscore_leg)
  
  # ── Node size legend ──────────────────────────────────────────────────────────
  legend_size_grob <- cowplot::get_legend(
    ggplot(data.frame(z = c(2, 4, 6)), aes(x = 1, y = z, size = z)) +
      geom_point(shape = 21, fill = "#aaaaaa", color = "black", stroke = 0.4) +
      scale_size_continuous(
        name   = "|Z|",
        breaks = c(2, 4, 6),
        range  = c(1.5, 5)
      ) +
      theme_void(base_size = 6) +
      theme(
        legend.position  = "left",
        legend.direction = "vertical",
        legend.title     = element_text(size = 6),
        legend.text      = element_text(size = 6),
        plot.margin      = margin(0, 0, 0, 0)
      )
  )
  
  # ── FDR significance note ─────────────────────────────────────────────────────
  legend_sig_grob <- ggplotGrob(
    ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
      geom_point(shape = 21, size = 4, fill = "grey70", color = "black") +
      geom_text(aes(label = "*"), color = "white", fontface = "bold", size = 3) +
      annotate("text", x = 1.3, y = 1, label = "FDR < 0.05",
               hjust = 0, size = 2) +
      xlim(0.8, 2) + ylim(0.5, 1.5) +
      theme_void(base_size = 6) +
      theme(plot.margin = margin(0, 0, 0, 0))
  )
  
  # ── Combine and save ──────────────────────────────────────────────────────────
  g_colorbars <- gridExtra::arrangeGrob(
    legend_scDRS_grob, legend_zscore_grob,
    ncol = 1, heights = c(1, 1)
  )
  legend_block <- gridExtra::arrangeGrob(
    g_colorbars, legend_size_grob, legend_sig_grob,
    ncol = 3, widths = c(1, 1, 1.2)
  )
  
  pdf(file   = file.path(resultdir, "WangNature/slingshot", cell_type,
                         paste0("lineage_scDRS_", version, "_",
                                gene_number, "genes_legend_", trait, ".pdf")),
      width  = 3.5, height = 4)
  grid::grid.newpage()
  grid::grid.draw(legend_block)
  dev.off()
  
  message("Saved legend: ", cell_type, " / ", version, " / ",
          gene_number, " genes / ", trait)
}


# ══════════════════════════════════════════════════════════════════════════════
# Preprocessing: compute scDRS mclust results and global ranges
# ══════════════════════════════════════════════════════════════════════════════

# ── whole gene set (run_v8, 253 genes) ────────────────────────────────────────
trait       <- "ASC_Pfdr"
version     <- "run_v8"
gene_number <- 253

# EN
ENlineage_scDRS_RVAS_sex_fdr001 <- fread(file.path(
  resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/EN",
  paste0(trait, ".scdrs_group.mclust")
)) %>%
  mutate(cluster = gsub("ENlineage_", "", group),
         cluster = as.character(cluster),
         label   = en_labels[cluster])

fwrite(ENlineage_scDRS_RVAS_sex_fdr001,
       file = file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/EN",
                        paste0(trait, ".scdrs_group.mclust.csv")),
       sep = "\t")

# IN
INlineage_scDRS_RVAS_sex_fdr001 <- fread(file.path(
  resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/IN",
  paste0(trait, ".scdrs_group.mclust")
)) %>%
  mutate(cluster = gsub("INlineage_", "", group),
         cluster = as.character(cluster),
         label   = in_labels[cluster])

fwrite(INlineage_scDRS_RVAS_sex_fdr001,
       file = file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/IN",
                        paste0(trait, ".scdrs_group.mclust.csv")),
       sep = "\t")

# Combine and compute FDR
lineage_scDRS_RVAS_sex_fdr001 <- rbind(
  ENlineage_scDRS_RVAS_sex_fdr001,
  INlineage_scDRS_RVAS_sex_fdr001
) %>%
  mutate(
    assoc_mcp_fdr  = p.adjust(assoc_mcp, method = "fdr"),
    assoc_mcp_bonf = p.adjust(assoc_mcp, method = "bonferroni")
  )

fwrite(lineage_scDRS_RVAS_sex_fdr001,
       file = file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes",
                        paste0(trait, ".scdrs_group.mclust.csv")),
       sep = "\t")

# Global z range
z_lim <- range(lineage_scDRS_RVAS_sex_fdr001$assoc_mcz, na.rm = TRUE)
saveRDS(z_lim, file.path(resultdir,
                         "WangNature/scDRS/run_v8/downstream_analysis/253_genes",
                         paste0("z_lim_", trait, ".rds")))

# ══════════════════════════════════════════════════════════════════════════════
# Load Seurat objects and slingshot curves (read once)
# ══════════════════════════════════════════════════════════════════════════════

dat_EN_lineage_mclust_filt <- readRDS(paste0(datapath, "dat_EN_lineage_mclust_slingshot.rds"))
dat_IN_lineage_mclust      <- readRDS(paste0(datapath, "dat_IN_lineage_mclust_slingshot.rds"))
crv_2d_EN <- readRDS(paste0(datapath, "curve_2d_EN_lineage.rds"))
crv_2d_IN <- readRDS(paste0(datapath, "curve_2d_IN_lineage.rds"))
sds_EN    <- SlingshotDataSet(crv_2d_EN)
sds_IN    <- SlingshotDataSet(crv_2d_IN)

# Load global ranges
z_lim <- readRDS(paste0(
  resultdir, "WangNature/scDRS/run_v8/downstream_analysis/253_genes/z_lim_", trait, ".rds"
))

# ══════════════════════════════════════════════════════════════════════════════
# Plot: whole gene set (run_v8, 253 genes)
# ══════════════════════════════════════════════════════════════════════════════
trait       <- "ASC_Pfdr"
version     <- "run_v8"
gene_number <- 253
lineage_scDRS_RVAS_sex_fdr001 <- fread(file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis",
                                                 paste0(gene_number,"_genes/",trait, ".scdrs_group.mclust.csv")))
p_EN <- plot_lineage(
  dat_seurat    = dat_EN_lineage_mclust_filt,
  sds           = sds_EN,
  lineage_scDRS = lineage_scDRS_RVAS_sex_fdr001,
  z_lim         = z_lim,
  trait         = trait,
  gene_number   = gene_number,
  version       = version,
  cell_type     = "EN",
  my_labels     = en_labels,
  resultdir     = resultdir
)

p_IN <- plot_lineage(
  dat_seurat    = dat_IN_lineage_mclust,
  sds           = sds_IN,
  lineage_scDRS = lineage_scDRS_RVAS_sex_fdr001,
  z_lim         = z_lim,
  trait         = trait,
  gene_number   = gene_number,
  version       = version,
  cell_type     = "IN",
  my_labels     = in_labels,
  resultdir     = resultdir
)

save_lineage_legend(
  p_scDRS     = p_EN,
  z_lim       = z_lim,
  version     = version,
  gene_number = gene_number,
  trait       = trait,
  resultdir   = resultdir,
  cell_type   = "EN"
)
