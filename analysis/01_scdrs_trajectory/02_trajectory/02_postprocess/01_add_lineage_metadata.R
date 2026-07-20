library(data.table)
library(R.utils)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(slingshot)
library(tidyverse)
library(scCustomize)


# ── Helper: add scDRS scores to Seurat object and generate plots ──────────────
#
# Reads a scDRS full_score file, attaches norm_score and nlog10_pval to the
# Seurat object as new metadata columns, and saves three figures:
#   1. UMAP FeaturePlot (PDF)
#   2. Lineage overlay on UMAP (PNG)
#   3. Smoothed norm_score vs pseudotime curve (PNG)
#
# @param seurat_obj   Seurat object to add metadata to
# @param version      scDRS run version, e.g. "run_v8" or "run_ddid_noddid"
# @param gene_number  gene set size or label, e.g. 253, 416, "ddid", "noddid"
# @param trait        trait identifier, e.g. "ASC_Pfdr"
# @param metadata_base  base metadata data.frame (cells × columns) with
#                       average_pseudotime and ID columns
# @param crv_2d       SlingshotDataSet object for lineage/curve overlay
# @param scoredir     base directory where scDRS score files live
# @param resultdir    base directory for saving output figures
# @param cell_type    "EN" or "IN" – determines subdirectory and reduction name
#
# @return the updated Seurat object (with new metadata columns added),
#         or the original object unchanged if the score file is missing
add_scDRS_and_plot <- function(seurat_obj, version, gene_number, trait,
                               metadata_base, crv_2d,
                               scoredir, resultdir,
                               cell_type = "EN") {
  
  # ── 1. Locate and load score file ───────────────────────────────────────────
  score_path <- file.path(
    scoredir, version, "score_file",
    paste0(gene_number, "_genes"),
    paste0(trait, ".full_score.gz")
  )
  
  if (!file.exists(score_path)) {
    message("  -> File not found, skipping: ", score_path)
    return(seurat_obj)   # return unchanged object
  }
  
  message("Processing ", version, " / ", gene_number, " genes / trait ", trait)
  
  score_df           <- as.data.frame(fread(score_path))
  rownames(score_df) <- score_df$V1
  
  scDRS_norm_score <- score_df %>%
    dplyr::select(norm_score, nlog10_pval)
  
  # Column names written into Seurat metadata
  score_col <- paste0("scDRS_", version, "_", gene_number, "_", trait, "_norm_score")
  pval_col  <- paste0("scDRS_", version, "_", gene_number, "_", trait, "_nlog10_pval")
  
  seurat_obj <- AddMetaData(
    seurat_obj,
    metadata = scDRS_norm_score,
    col.name = c(score_col, pval_col)
  )
  
  # ── 2. Shared output prefix ─────────────────────────────────────────────────
  out_dir    <- file.path(resultdir, "WangNature/slingshot", cell_type)
  # Reduction name follows the pattern wnn.<lower cell_type>.umap2
  reduction  <- paste0("wnn.", cell_type, ".umap2")
  file_stem  <- paste0("scDRS_", version, "_", gene_number, "genes_", trait)
  
  # ── 3. UMAP FeaturePlot (PDF) ────────────────────────────────────────────────
  p_scDRS <- FeaturePlot_scCustom(
    seurat_obj,
    features   = score_col,
    reduction  = reduction,
    colors_use = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                   "#F7F7F7",
                   "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    max.cutoff = 5,
    min.cutoff = -3,
    na_cutoff  = NULL,
    raster     = TRUE
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stem, "_norm_score_umap.pdf")),
    plot     = p_scDRS,
    width = 7, height = 5, units = "in"
  )
  
  # Extract point coordinates and colours for base-R overlay plots
  p_build  <- ggplot_build(p_scDRS)
  col_vec  <- unlist(p_build$data[[1]]["colour"])
  pts      <- c(p_build$data[[1]]["x"], p_build$data[[1]]["y"])
  umap_labs <- c(paste0(reduction, "_1"), paste0(reduction, "_2"))
  
  # ── Internal helper: draw scatter + slingshot overlay, save as PNG ──────────
  save_slingshot_png <- function(filename, sling_type) {
    png(filename = filename, width = 8, height = 6, units = "in", res = 300)
    plot(pts,
         xlab = umap_labs[1], ylab = umap_labs[2],
         col  = col_vec, asp = 1, pch = 16, cex = 0.2)
    lines(SlingshotDataSet(crv_2d),
          lwd = 1.5, col = "black",
          type = sling_type,
          show.constraints = TRUE)
    dev.off()
  }
  
  # ── 4. Lineage overlay (PNG) ─────────────────────────────────────────────────
  save_slingshot_png(
    filename   = file.path(out_dir, paste0("lineage_", file_stem, "_norm_score_umap.png")),
    sling_type = "lineages"
  )
  
  # ── 5. Curve overlay (PNG) ───────────────────────────────────────────────────
  save_slingshot_png(
    filename   = file.path(out_dir, paste0("curve_", file_stem, "_norm_score_umap.png")),
    sling_type = "curve"
  )
  
  # ── 6. Smoothed norm_score vs pseudotime (PNG) ───────────────────────────────
  scDRS_smooth <- metadata_base %>%
    left_join(
      scDRS_norm_score %>% mutate(ID = rownames(scDRS_norm_score)),
      by = "ID"
    ) %>%
    dplyr::select(ID, norm_score, average_pseudotime) %>%
    mutate(bin = cut(average_pseudotime, breaks = 200)) %>%
    group_by(bin) %>%
    summarise(
      average_pseudotime = mean(average_pseudotime),
      scDRS_norm_score   = mean(norm_score),
      .groups = "drop"
    )
  
  p_curve <- ggplot(scDRS_smooth,
                    aes(x = average_pseudotime, y = scDRS_norm_score)) +
    geom_line(color = "#6a3d9a", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_classic()
  
  ggsave(
    filename = file.path(out_dir, paste0(file_stem, "_norm_score_avgpseudotime.png")),
    plot     = p_curve,
    width = 7, height = 5, units = "in"
  )
  
  invisible(seurat_obj)   # return updated object
}


# ── Main loop ─────────────────────────────────────────────────────────────────

# EN ----
datapath  <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
scoredir  <- file.path(resultdir, "WangNature/scDRS")

dat_EN_lineage_mclust_filt <- readRDS(paste0(datapath, "dat_EN_lineage_mclust_slingshot.rds"))
crv_2d_EN                  <- readRDS(paste0(datapath, "curve_2d_EN_lineage.rds"))
metadata_EN              <- dat_EN_lineage_mclust_filt@meta.data

gene_numbers <- c(253, 416, 696, 951, "ddid", "noddid",285, 100, 200)
versions     <- c("run_v8", "run_ddid_noddid", "run_GWAS_v3")
traits       <- c("ASC_Pfdr",
                  "ASC_Pfdr_ddid_new_DMN", "ASC_Pfdr_noddid_new_DMN", "ASC_Pfdr_ddid_perproband_new_DMN", "ASC_Pfdr_noddid_perproband_new_DMN", "PGC_ASD_2019", "ASD_Matoba_2020")

# Iterate over all combinations; seurat object is updated in-place each round
for (version in versions) {
  for (gene_number in gene_numbers) {
    for (trait in traits) {
      dat_EN_lineage_mclust_filt <- add_scDRS_and_plot(
        seurat_obj    = dat_EN_lineage_mclust_filt,
        version       = version,
        gene_number   = gene_number,
        trait         = trait,
        metadata_base = metadata_EN,
        crv_2d        = crv_2d_EN,
        scoredir      = scoredir,
        resultdir     = resultdir,
        cell_type     = "EN"
      )
    }
  }
}


# IN ----
datapath  <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
scoredir  <- file.path(resultdir, "WangNature/scDRS")
dat_IN_lineage_mclust <- readRDS(paste0(datapath, "dat_IN_lineage_mclust_slingshot.rds"))
crv_2d_IN                  <- readRDS(paste0(datapath, "curve_2d_IN_lineage.rds"))
metadata_base_IN           <- dat_IN_lineage_mclust@meta.data

gene_numbers <- c(253, 416, 696, 951, "ddid", "noddid",285, 100, 200)
versions     <- c("run_v8", "run_ddid_noddid", "run_GWAS_v3")
traits       <- c("ASC_Pfdr",
                  "ASC_Pfdr_ddid_new_DMN", "ASC_Pfdr_noddid_new_DMN", "ASC_Pfdr_ddid_perproband_new_DMN", "ASC_Pfdr_noddid_perproband_new_DMN", "PGC_ASD_2019", "ASD_Matoba_2020")

for (version in versions) {
  for (gene_number in gene_numbers) {
    for (trait in traits) {
      dat_IN_lineage_mclust <- add_scDRS_and_plot(
        seurat_obj    = dat_IN_lineage_mclust,
        version       = version,
        gene_number   = gene_number,
        trait         = trait,
        metadata_base = metadata_base_IN,
        crv_2d        = crv_2d_IN,
        scoredir      = scoredir,
        resultdir     = resultdir,
        cell_type     = "IN"        
      )
    }
  }
}

# save ----
saveRDS(dat_EN_lineage_mclust_filt, paste0(datapath, "dat_EN_lineage_mclust_slingshot.rds"))
metadata_final <- data.frame(dat_EN_lineage_mclust_filt[[]])
write.csv(metadata_final, file = paste0(resultdir,"WangNature/slingshot/EN/EN_lineage_slingshot_results.csv"),row.names = TRUE)

saveRDS(dat_IN_lineage_mclust, paste0(datapath, "dat_IN_lineage_mclust_slingshot.rds"))
metadata_final_IN <- data.frame(dat_IN_lineage_mclust[[]])
write.csv(metadata_final_IN, file = paste0(resultdir,"WangNature/slingshot/IN/IN_lineage_slingshot_results.csv"),row.names = TRUE)

