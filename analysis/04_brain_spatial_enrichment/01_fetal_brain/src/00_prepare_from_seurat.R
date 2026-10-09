#!/usr/bin/env Rscript
# =============================================================================
# Fig2-occipital-only: stage 0 of the pipeline (src/00-08)
#
# Prepares the GW20 occipital Visium section for scDRS scoring. Exports the raw
# count matrix, the spot coordinates and the authors' own laminar and areal
# annotation, which later stages use directly - there is no label transfer and
# no integration step anywhere in this analysis.
#
# Input    one deposited Seurat object from Zenodo 10.5281/zenodo.14422018
#            Visium_A1_brain_011124.rds   GW20 occipital (V1+V2), 3,591 spots  309,512,278 B
# Outputs  export/A1_counts.mtx + features/barcodes   -> scDRS input
#          export/A1_author_labels.csv                -> authors' per-spot annotation
#          export/A1_coords_scales.parquet            -> spot coordinates
#          export/A1_histology_lowres.png             -> H&E image
#
# A prefrontal section (D1_brain_no15.rds) from the same deposit was scored
# alongside this one in an earlier version of the analysis. It was dropped: the
# two sections come from different donors, so no prefrontal-versus-occipital
# contrast is separable from a donor difference. Nothing in this pipeline reads
# it, and the joint CCA/Harmony clustering that transferred layer labels to it
# was removed with it.
#
# R packages used: Seurat 5.5.1, SeuratObject 5.4.0, Matrix 1.7.6, arrow 25.0.0,
#                  dplyr 1.2.1, png 0.1.9 (full list: environment_r.txt)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(arrow)
  library(dplyr); library(png)
})
set.seed(1)

# Directory holding the two deposited Seurat objects (Zenodo
# 10.5281/zenodo.14422018). Set VISIUM_RDS_DIR to point at it; the default is
# the working directory, so the script runs unedited from there.
ROOT <- Sys.getenv("VISIUM_RDS_DIR", unset = ".")
OBJ  <- list(A1 = "Visium_A1_brain_011124.rds")

absent_files <- Filter(function(f) !file.exists(file.path(ROOT, f)), unlist(OBJ))
if (length(absent_files)) {
  stop("Deposited Visium objects not found in '", ROOT, "': ",
       paste(absent_files, collapse = ", "),
       "\nDownload them from Zenodo 10.5281/zenodo.14422018 and set",
       " VISIUM_RDS_DIR to that directory.")
}
# name of each section's image within its Seurat object
IMG  <- list(A1 = "slice1")
dir.create("export", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. The authors' laminar and areal annotation of the occipital section
#
# The deposited object does not carry it (meta.data has 8 columns, neither
# assay has feature metadata, and @misc and @tools are empty). The annotation
# exists only as a RenameIdents call in Fig4/Visium_A1_Fig4_h_i_j.R of the
# authors' analysis repository, which their Code Availability statement names:
#   github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-
#   Layer-and-Area-Specification
# AUTHOR_MAP below is transcribed from that call and matches the repository at
# commit 0f5fa04451d2ccac63bb01e26184cfc65162cfff. The call renames the clusters
# of seurat_clusters (= SCT_snn_res.1.5, levels 0-18). Each new name gives the
# layer and, where the annotation is area-specific, the cortical area (V1 or V2).
# Cluster 14 (128 spots) has no entry: the authors dropped it, and their
# downstream script ED Fig18/EDFig18D_sc_visium.R reads
# "Visium_A1brain_no14_annotated.rds".
# The prefrontal section has no equivalent annotation; its figure script,
# ED Fig5/Visium_D1_EDFig5_g_h_i_j.R, contains no RenameIdents.
# -----------------------------------------------------------------------------
AUTHOR_MAP <- c(
  "0" = "iSVZ",   "1" = "IZ",      "2" = "oSVZ-1",  "3" = "VZ",     "4" = "oSVZ-2",
  "5" = "SP-V1",  "6" = "oSVZ-L1", "7" = "oSVZ-L2", "8" = "L4-V2",  "9" = "L5/6-V1",
  "10" = "L4-V1", "11" = "oSVZ-3", "12" = "L2/3-V2","13" = "SP-V2", "15" = "L5/6-V2",
  "16" = "L2-V1", "17" = "L3-V1",  "18" = "iSVZ-L")

# -----------------------------------------------------------------------------
# 2. Per-section export: raw counts for scDRS, spot coordinates, histology image
#
# Spot coordinates. In these objects imagerow/imagecol from GetTissueCoordinates
# are already in pixels of the low-resolution (600x600) image, so they must not
# be multiplied by scale.factors$lowres (0.2). Checked by sampling eosin
# intensity at the spot positions: with the coordinates as returned, 100.0% of
# spots fall on tissue; scaled by 0.2, only 44.7% do.
# Orientation. imagerow increases toward the pia, so plotting y = imagerow on a
# non-inverted axis gives the published layout, cortical plate above and
# germinal zones below. The authors' own script flips the image instead, with
# rotateSeuratImage(rotation = "180").
# -----------------------------------------------------------------------------
for (tag in names(OBJ)) {
  visium_obj <- readRDS(file.path(ROOT, OBJ[[tag]]))
  DefaultAssay(visium_obj) <- "Spatial"

  # raw counts, not SCT-normalised values: scDRS normalises them itself
  section_counts <- LayerData(visium_obj, assay = "Spatial", layer = "counts")
  Matrix::writeMM(section_counts, sprintf("export/%s_counts.mtx", tag))
  write.table(rownames(section_counts), sprintf("export/%s_features.tsv", tag),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(section_counts), sprintf("export/%s_barcodes.tsv", tag),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  im  <- visium_obj@images[[IMG[[tag]]]]
  tissue_coords <- GetTissueCoordinates(visium_obj, image = IMG[[tag]])
  coord_table <- data.frame(barcode  = rownames(tissue_coords),
                            row_raw  = tissue_coords$imagerow,   # already in low-res pixels
                            col_raw  = tissue_coords$imagecol,
                            # the two scaled columns are kept only to document the
                            # rejected convention
                            row_lowres = tissue_coords$imagerow * im@scale.factors$lowres,
                            col_lowres = tissue_coords$imagecol * im@scale.factors$lowres)
  write_parquet(coord_table, sprintf("export/%s_coords_scales.parquet", tag))
  png::writePNG(im@image, sprintf("export/%s_histology_lowres.png", tag))

  if (tag == "A1") {
    stopifnot(identical(levels(Idents(visium_obj)), as.character(0:18)))
    spot_label <- as.character(AUTHOR_MAP[as.character(visium_obj$seurat_clusters)])
    spot_label[is.na(spot_label)] <- "14"        # the authors' dropped cluster
    label_table <- data.frame(barcode = colnames(visium_obj),
                              cluster = as.character(visium_obj$seurat_clusters),
                              author_label = spot_label,
                              imagerow = tissue_coords[colnames(visium_obj), "imagerow"],
                              imagecol = tissue_coords[colnames(visium_obj), "imagecol"],
                              dropped_by_authors = spot_label == "14")
    write.csv(label_table, "export/A1_author_labels.csv", row.names = FALSE)
    stopifnot(sum(label_table$dropped_by_authors) == 128, nrow(label_table) == 3591)
  }
  rm(visium_obj, section_counts, im); invisible(gc())
}

# -----------------------------------------------------------------------------
cat("stage 0 complete\n")
