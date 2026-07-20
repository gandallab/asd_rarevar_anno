library(Seurat)
library(Signac)
library(slingshot)
library(tidyverse)
library(scCustomize)
library(mclust)
library(future)
library(patchwork)
library(dplyr)

plan("multicore", workers = 12)
# options(future.globals.maxSize = 2000000 * 1024^2)

datapath = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
obj = readRDS(paste0(datapath, "snMultiome_atlas_Seurat_object.rds")) # 74GB

# IPC.Glia - Tri-IPC to Astro/OPC/IN
# focus on IN lineage -----
dat_IN_lineage <- subset(obj, subset = type %in% c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-Glia",
                                                   "IN-dLGE-Immature","IN-CGE-Immature", 
                                                   "IN-CGE-VIP", "IN-CGE-SNCG", 
                                                   "IN-Mix-LAMP5","IN-MGE-Immature", 
                                                   "IN-MGE-SST", "IN-MGE-PV"))

# release memory
rm(obj)

#WNN and UMAP
dat_IN_lineage <- FindMultiModalNeighbors(dat_IN_lineage,
                                          reduction.list = list("pca", "integrated_lsi"),
                                          dims.list = list(1:50, 2:40),
                                          knn.graph.name = "wknn.IN",
                                          snn.graph.name = "wsnn.IN",
                                          weighted.nn.name = "weighted.nn.IN")

dat_IN_lineage <- RunUMAP(dat_IN_lineage,
                          nn.name = "weighted.nn.IN",
                          reduction.name = "wnn.IN.umap",
                          reduction.key = "wnninumap_",
                          min.dist = 0.3,
                          n.neighbors = 50L,
                          n.components = 8)
dat_IN_lineage <- RunUMAP(dat_IN_lineage,
                          nn.name = "weighted.nn.IN",
                          reduction.name = "wnn.IN.umap2",
                          reduction.key = "wnninumap2_",
                          min.dist = 0.3,
                          n.neighbors = 50L,
                          n.components = 2)

#save data
saveRDS(dat_IN_lineage, paste0(datapath, "dat_IN_lineage.rds"))


#plot metadata
DimPlot_scCustom(dat_IN_lineage,
                 group.by = "subclass",
                 reduction = "wnn.IN.umap2",
                 raster = FALSE)
ggsave(file = paste0(resultdir,"WangNature/slingshot/IN/subclass_umap.png"),width = 8, height = 5, units = "in")

DimPlot_scCustom(dat_IN_lineage,
                 group.by = "type",
                 reduction = "wnn.IN.umap2",
                 colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C"), raster = FALSE)
ggsave(file = paste0(resultdir,"WangNature/slingshot/IN/type_umap.png"),width = 8, height = 5, units = "in")

DimPlot_scCustom(dat_IN_lineage,
                 group.by = "dataset",
                 reduction = "wnn.IN.umap2",
                 raster = FALSE)
ggsave(file = paste0(resultdir,"WangNature//slingshot/IN/dataset_umap_filtered.png"),width = 15, height = 9, units = "in")

DimPlot_scCustom(dat_IN_lineage,
                 group.by = "region_summary",
                 reduction = "wnn.IN.umap2",
                 colors_use = c("#009e73", "#ffa500", "#0072b2"), raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/region_summary.png"),width = 7, height = 5, units = "in")

DimPlot_scCustom(dat_IN_lineage,
                 group.by = "Group",
                 reduction = "wnn.IN.umap2",
                 colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/group_umap_filtered.png"),width = 7, height = 5, units = "in")

FeaturePlot_scCustom(dat_IN_lineage,
                     features = "log2_age",
                     reduction = "wnn.IN.umap2",
                     na_cutoff = NULL,
                     raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/log2_age_umap_filtered.png"),width = 7, height = 5, units = "in")

# get different levels of clusters using mclust clusters -----
rd = Embeddings(dat_IN_lineage, reduction="wnn.IN.umap")
dat_IN_lineage_mclust <- dat_IN_lineage
rm(dat_IN_lineage)

for (i in c(14:23)) {
  set.seed(2025)
  cl <- Mclust(rd, G = i)$classification
  dat_IN_lineage_mclust <-  AddMetaData(dat_IN_lineage_mclust,
                                        metadata = cl,
                                        col.name = paste0("mclust",i))
  DimPlot_scCustom(dat_IN_lineage_mclust,
                   group.by = paste0("mclust",i),
                   reduction = "wnn.IN.umap2",
                   label = TRUE,
                   raster = FALSE)
  ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/mclust/mclust",i, "_umap.png"), width = 8, height = 5, units = "in")
}

#save data with mclust information
saveRDS(dat_IN_lineage_mclust, paste0(datapath, "dat_IN_lineage_mclust.rds"))

# perform slingshot trajectory inference-----
rd = Embeddings(dat_IN_lineage_mclust, reduction = "wnn.IN.umap")
cl <- dat_IN_lineage_mclust$mclust22
times <- dat_IN_lineage_mclust$log2_age

# merge two clusters belong to same cell type
cl_merged <- cl
cl_merged[cl == "21"] <- "15"

lin <- getLineages(
  rd, cl_merged,
  start.clus = "19", #RG-vRG
  end.clus = c("16","20","5","17","10","18"),
  times = times,
  use.median = FALSE,
  dist.method = "slingshot"
)
slingLineages(lin)

crv <- getCurves(lin, approx_points = 200)
saveRDS(crv, paste0(datapath, "curve_IN_lineage.rds"))

pseudotime_mtx <- slingPseudotime(crv) # pseudotime values for each lineage (cell × lineages matrix)
pseudotime_avg <- slingAvgPseudotime(crv) # average pseudotime per cell (across all relevant lineages)
cell_weights <- slingCurveWeights(crv) # lineage weights per cell (probability of belonging to each lineage)
branchID <- slingBranchID(crv) # branch assignment for each cell

rd_2d <- Embeddings(dat_IN_lineage_mclust, reduction = "wnn.IN.umap2")
crv_2d <- embedCurves(crv, rd_2d) # project the slingshot curves into 2D space
saveRDS(crv_2d, paste0(datapath, "curve_2d_IN_lineage.rds"))

dat_IN_lineage_mclust$average_pseudotime <- pseudotime_avg
colnames(pseudotime_mtx) <- c(
  "IN.MGE.PV.2_pseudotime",
  "IN.MGE.SST.1_pseudotime",
  "IN.MGE.PV.1_pseudotime",
  "IN.CGE.VIP_pseudotime",
  "IN.MGE.SST.2_pseudotime",
  "IN.Mix.LAMP5_pseudotime",
  "IN.CGE.SNCG_pseudotime",
  "RG_pseudotime",
  "IPC.Glia_pseudotime"
)
colnames(cell_weights) <- c(
  "IN.MGE.PV.2_lineageWeight",
  "IN.MGE.SST.1_lineageWeight",
  "IN.MGE.PV.1_lineageWeight",
  "IN.CGE.VIP_lineageWeight",
  "IN.MGE.SST.2_lineageWeight",
  "IN.Mix.LAMP5_lineageWeight",
  "IN.CGE.SNCG_lineageWeight",
  "RG_lineageWeight",
  "IPC.Glia_lineageWeight"
)

dat_IN_lineage_mclust <- AddMetaData(
  dat_IN_lineage_mclust,
  metadata = pseudotime_mtx,
  col.name = colnames(pseudotime_mtx)
)
dat_IN_lineage_mclust <- AddMetaData(
  dat_IN_lineage_mclust,
  metadata = cell_weights,
  col.name = colnames(cell_weights)
)

dat_IN_lineage_mclust@meta.data = dat_IN_lineage_mclust@meta.data %>% mutate(mclust = paste0("INlineage_", mclust22))
metadata_mclust <- data.frame(dat_IN_lineage_mclust[[]]) %>% select("ID",	"type",	"mclust22",	"mclust")
write.csv(metadata_mclust, file = paste0(resultdir,"WangNature/slingshot/IN/IN_lineage_mclust.csv"))
metadata_final <- data.frame(dat_IN_lineage_mclust[[]])
write.csv(metadata_final, file = paste0(resultdir,"WangNature/slingshot/IN/IN_lineage_slingshot_results.csv"),row.names = TRUE)

saveRDS(dat_IN_lineage_mclust, paste0(datapath, "dat_IN_lineage_mclust_slingshot.rds"))


#plot slingshot results -----
DimPlot_scCustom(dat_IN_lineage_mclust,
                 group.by = "mclust22",
                 colors_use = DiscretePalette_scCustomize(num_colors = 22,
                                                          palette = "alphabet2"),
                 reduction = "wnn.IN.umap2",
                 label = TRUE,
                 raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/mclust/mclust22_umap.png"), width = 8, height = 5, units = "in")

#lineage with mclust
png(filename=paste0(resultdir,"WangNature/slingshot/IN/lineage_mclust22_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(rd_2d,
     col = DiscretePalette_scCustomize(num_colors = 21,
                                       palette = "alphabet2")[cl_merged],
     asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()

#curve with mclust
png(filename=paste0(resultdir,"WangNature/slingshot/IN/curve_mclust22_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(rd_2d,
     col = DiscretePalette_scCustomize(num_colors = 21, palette = "alphabet2")[cl_merged],
     asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with type
p_type <-  DimPlot_scCustom(dat_IN_lineage_mclust,
                            group.by = "type",
                            colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
                                           "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C"),
                            reduction = "wnn.IN.umap2",
                            raster = FALSE)
p_type <- ggplot_build(p_type)
col <- unlist(p_type$data[[1]]["colour"])

png(filename=paste0(resultdir,"WangNature/slingshot/IN/lineage_type_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_type$data[[1]]["x"],
       p_type$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()

#curve with type
png(filename=paste0(resultdir,"WangNature/slingshot/IN/curve_type_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_type$data[[1]]["x"],
       p_type$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#curve with type more transparent
png(filename=paste0(resultdir,"WangNature/slingshot/IN/curve_type_umap_transparent.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_type$data[[1]]["x"],
       p_type$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = scales::alpha(col,0.1), asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with age group
p_age_group <-  DimPlot_scCustom(dat_IN_lineage_mclust,
                                 group.by = "Group",
                                 colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"),
                                 reduction = "wnn.IN.umap2",
                                 raster = FALSE)
p_age_group <- ggplot_build(p_age_group)
col <- unlist(p_age_group$data[[1]]["colour"])

png(filename=paste0(resultdir,"WangNature/slingshot/IN/lineage_age_group_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_age_group$data[[1]]["x"],
       p_age_group$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()

#curve with age group
png(filename=paste0(resultdir,"WangNature/slingshot/IN/curve_age_group_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_age_group$data[[1]]["x"],
       p_age_group$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with region
p_region <-  DimPlot_scCustom(dat_IN_lineage_mclust,
                              group.by = "region_summary",
                              colors_use = c("#009e73", "#ffa500", "#0072b2"),
                              reduction = "wnn.IN.umap2",
                              raster = FALSE)
p_region <- ggplot_build(p_region)
col <- unlist(p_region$data[[1]]["colour"])
png(filename=paste0(resultdir,"WangNature/slingshot/IN/lineage_region_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_region$data[[1]]["x"],
       p_region$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()

#curve with age group
png(filename=paste0(resultdir,"WangNature/slingshot/IN/curve_region_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_region$data[[1]]["x"],
       p_region$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 1.5, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with pseudotime
p_psueodtime <- FeaturePlot_scCustom(dat_IN_lineage_mclust,
                                     features = "average_pseudotime",
                                     reduction = "wnn.IN.umap2",
                                     colors_use = viridis_light_high,
                                     na_cutoff = NULL,
                                     raster = T)
ggsave(file =paste0(resultdir,"WangNature/slingshot/IN/average_pseudotime.pdf"),
       width = 7, height = 5, units = "in")

p_psueodtime <- ggplot_build(p_psueodtime)
col <- unlist(p_psueodtime$data[[1]]["colour"])

png(filename=paste0(resultdir,"WangNature/slingshot/IN/lineage_pseudotime_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_psueodtime$data[[1]]["x"],
       p_psueodtime$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2",
     col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()

#curve with pseudotime
png(file=paste0(resultdir,"WangNature/slingshot/IN/curve_pseudotime_umap.png"), width = 8, height = 6, units = 'in', res = 300)
plot(c(p_psueodtime$data[[1]]["x"],
       p_psueodtime$data[[1]]["y"]),
     xlab="wnninumap2_1",
     ylab="wnninumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d),
      lwd = 1.5, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#plot pseudotime in individual lineages
for (i in c("IN.MGE.PV.2_pseudotime",
            "IN.MGE.SST.1_pseudotime",
            "IN.MGE.PV.1_pseudotime",
            "IN.CGE.VIP_pseudotime",
            "IN.MGE.SST.2_pseudotime",
            "IN.Mix.LAMP5_pseudotime",
            "IN.CGE.SNCG_pseudotime",
            "RG_pseudotime",
            "IPC.Glia_pseudotime")) {
  p <- FeaturePlot_scCustom(dat_IN_lineage_mclust,
                            features = i,
                            reduction = "wnn.IN.umap2",
                            colors_use = viridis_light_high,
                            na_cutoff = NULL, raster = FALSE)
  ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/individual_lineages/", i, ".png"),
         width = 7, height = 5, units = "in")
}

#psuedotime of all lineages in the same plot
pseudotime_list <- list()
for (i in c("IN.MGE.PV.2_pseudotime",
            "IN.MGE.SST.1_pseudotime",
            "IN.MGE.PV.1_pseudotime",
            "IN.CGE.VIP_pseudotime",
            "IN.MGE.SST.2_pseudotime",
            "IN.Mix.LAMP5_pseudotime",
            "IN.CGE.SNCG_pseudotime",
            "RG_pseudotime",
            "IPC.Glia_pseudotime")) {
  p <- FeaturePlot_scCustom(dat_IN_lineage_mclust,
                            features = i,
                            reduction = "wnn.IN.umap2",
                            colors_use = viridis_light_high,
                            na_cutoff = NULL, raster = FALSE) +
    labs(title = NULL) + NoAxes() + NoLegend() +
    scale_colour_viridis_c(limits = c(0,26.48495),
                           oob=scales::squish, na.value = "lightgrey")
  pseudotime_list[[i]] <- p
}
patchwork::wrap_plots(pseudotime_list, nrow = 2)
ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/individual_lineages/all_pseudotime.png"),width = 20, height = 7, units = "in")

#highlight cells in individual lineages
colors <- c('#00429d', '#3e67ae', '#618fbf', '#85b7ce', '#b1dfdb', '#ffcab9', '#fd9291', '#e75d6f', '#8bc53f')
lineages <- c("IN.MGE.PV.2", "IN.MGE.SST.1", "IN.MGE.PV.1", "IN.CGE.VIP", "IN.MGE.SST.2", "IN.Mix.LAMP5", "IN.CGE.SNCG", "RG", "IPC.Glia")
for (i in 1:9) {
  Idents(dat_IN_lineage_mclust) <- ifelse(is.na(dat_IN_lineage_mclust[[paste0(lineages[i], "_pseudotime")]]), "ident.remove", "ident.keep")
  cells <- WhichCells(object = dat_IN_lineage_mclust, ident = "ident.keep")
  cells <- list(cells)
  names(cells) <- lineages[i]
  p <- Cell_Highlight_Plot(dat_IN_lineage_mclust, cells_highlight = cells, reduction = "wnn.IN.umap2", highlight_color = colors[i], raster = FALSE)
  assign(lineages[i], p)
  ggsave(file=paste0(resultdir,"WangNature/slingshot/IN/individual_lineages/", lineages[i], "_assignment.png"),width = 7, height = 5, units = "in")
}

#save plot regarding "IN.MGE.PV.1","IN.MGE.PV.2", "IN.MGE.SST.1", "IN.MGE.SST.2", "IN.CGE.VIP", "IN.Mix.LAMP5", "IN.CGE.SNCG"
IN.MGE.PV.1 + IN.MGE.PV.2 + IN.MGE.SST.1 + IN.MGE.SST.2 + IN.CGE.VIP + IN.Mix.LAMP5 + IN.CGE.SNCG + plot_layout(ncol = 2)
ggsave(file = paste0(resultdir,"WangNature/slingshot/IN/individual_lineages/BP2.png"), width = 23, height = 8, units = "in")
