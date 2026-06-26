library(Seurat)
library(SeuratDisk)
library(dplyr)
library(data.table)
library(ggplot2)
library(readxl)
datapath = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
datadir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/"

# A Seurat object containing count matrices, dimension reduction, and metadata  from 38 samples in this experiment (snMultiome_atlas_Seurat_object.rds) was downloaded from https://datadryad.org/dataset/doi:10.5061/dryad.2280gb612#citations.
obj = readRDS(paste0(datapath, "snMultiome_atlas_Seurat_object.rds")) # 74GB

DefaultAssay(obj) <- "RNA"
obj_rna_raw <- DietSeurat(
  obj,
  assays    = c("RNA","SCT"),
  dimreducs = NULL,
  graphs    = NULL
)
DefaultAssay(obj_rna_raw) <- "RNA"

# performed size factor normalization (10,000 counts per cell) and log transformation (log1p)
# As I've already prepared the snRNA data, when compute scDRS score and perform downstream analysis, use --flag-raw-count False - if didn't do normalization here, use --flag-raw-count True
obj_rna_raw <- NormalizeData(
  object = obj_rna_raw,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
obj_rna_raw$RNA@data

# transfer factor col in meta.data to character (including type, avioding 0/1/2 code in AnnData)
f_is_factor <- sapply(obj_rna_raw@meta.data, is.factor)
obj_rna_raw@meta.data[f_is_factor] <- lapply(
  obj_rna_raw@meta.data[f_is_factor],
  as.character
)
obj_rna_raw@meta.data$region_combine = obj_rna_raw@meta.data$Region
obj_rna_raw@meta.data$region_combine[obj_rna_raw@meta.data$region_combine == "V1"] = "BA17"

# transfer to h5Seurat
SaveH5Seurat(
  obj_rna_raw,
  filename = paste0(datapath, "obj_rna_raw.h5seurat"),
  overwrite = TRUE
)

# transfer to h5ad
Convert(
  paste0(datapath, "obj_rna_raw.h5seurat"),
  dest      = "h5ad",
  assay     = "RNA",     # same as DefaultAssay
  overwrite = TRUE
)