library(Seurat)
library(Signac)
library(future)
library(tidyverse)
library(scater)
library(scran)
library(edgeR)
library(glmmSeq)
library(scCustomize)
library(patchwork)
library(ggrepel)
library(data.table)
library(dplyr)
library(anndata)
library(reticulate)
# set up virtual env every time
use_virtualenv("py310-anndata", required = TRUE)
py_config()

# plan("multisession", workers = 12)
# options(future.globals.maxSize = 2000000 * 1024^2)

#import primary data
datapath = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
trait = "ASC_Pfdr" 
dat_IN_lineage_mclust <- readRDS(paste0(datapath, "dat_IN_lineage_mclust_slingshot.rds"))

#subset to focus on MGE_SST
primary_MGE_SST <- subset(dat_IN_lineage_mclust, subset = type == "IN-MGE-SST")
DefaultAssay(primary_MGE_SST) <- "integrated"

#re-annotate primary data
primary_MGE_SST <- RunUMAP(primary_MGE_SST, reduction = "pca", dims = 1:50)
primary_MGE_SST <- FindNeighbors(primary_MGE_SST, dims = 1:50)
primary_MGE_SST <- FindClusters(primary_MGE_SST, resolution = 0.5)
DimPlot_scCustom(primary_MGE_SST, reduction = "umap", group.by = "integrated_snn_res.0.5", raster = FALSE, label = T)
ggsave(filename = paste0(resultdir,
                         "WangNature/MGE_SST/primary_MGE_SST_integrated_snn_res.0.5.pdf"), width = 8, height = 6)

DimPlot_scCustom(primary_MGE_SST, reduction = "umap", group.by = "mclust22", raster = FALSE, label = T)
ggsave(filename = paste0(resultdir,
                         "WangNature/MGE_SST/primary_MGE_SST_mclust22.pdf"), width = 8, height = 6)

# subtype
primary_MGE_SST@meta.data$subtype = NA
primary_MGE_SST@meta.data$subtype[primary_MGE_SST@meta.data$mclust22==3 | primary_MGE_SST@meta.data$mclust22==5] = "IN-MGE-SST-1"
primary_MGE_SST@meta.data$subtype[primary_MGE_SST@meta.data$mclust22==9 | primary_MGE_SST@meta.data$mclust22==11] = "IN-MGE-SST-2"
primary_MGE_SST@meta.data$subtype[is.na(primary_MGE_SST@meta.data$subtype)] = "Other"
DimPlot_scCustom(primary_MGE_SST, reduction = "umap", group.by = "subtype", raster = FALSE, label = T, colors_use = c("#e63946", "#1d3557","grey"))
ggsave(filename = paste0(resultdir,
                         "WangNature/MGE_SST/primary_MGE_SST_subtype.pdf"), width = 8, height = 6)

#save data
saveRDS(primary_MGE_SST, paste0(datapath, "primary_MGE_SST.rds"))
#plot metadata
DimPlot_scCustom(primary_MGE_SST,
                 group.by = "subclass",
                 reduction = "umap",
                 raster = FALSE)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/subclass_umap.png"),width = 8, height = 5, units = "in")

DimPlot_scCustom(primary_MGE_SST,
                 group.by = "type",
                 reduction = "umap",
                 colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C"), raster = FALSE)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/type_umap.png"),width = 8, height = 5, units = "in")

DimPlot_scCustom(primary_MGE_SST,
                 group.by = "dataset",
                 reduction = "umap",
                 raster = FALSE)
ggsave(file = paste0(resultdir,"WangNature//slingshot/IN/dataset_umap_filtered.png"),width = 15, height = 9, units = "in")

DimPlot_scCustom(primary_MGE_SST,
                 group.by = "region_summary",
                 reduction = "umap",
                 colors_use = c("#009e73", "#ffa500", "#0072b2"), raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/MGE_SST/region_summary.pdf"),width = 7, height = 5, units = "in")

DimPlot_scCustom(primary_MGE_SST,
                 group.by = "Group",
                 reduction = "umap",
                 colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/MGE_SST/group_umap_filtered.pdf"),width = 7, height = 5, units = "in")

FeaturePlot_scCustom(primary_MGE_SST,
                     features = "log2_age",
                     reduction = "umap",
                     na_cutoff = NULL,
                     raster = FALSE)
ggsave(file=paste0(resultdir,"WangNature/MGE_SST/log2_age_umap_filtered.pdf"),width = 7, height = 5, units = "in")

# MapMyCell -----
dataIn <- primary_MGE_SST@assays$SCT$counts 
annoIn <- primary_MGE_SST@meta.data
dataQC <- dataIn[,colSums(dataIn)>250] # all barcodes with >250 reads as interesting "cells"
dim(dataIn)
dim(dataQC)

# output to h5ad format
dataQCt = Matrix::t(dataQC)
# Convert to anndata format
ad <- AnnData(
  X = dataQCt,
  obs = data.frame(group = rownames(dataQCt), row.names = rownames(dataQCt)),
  var = data.frame(type = colnames(dataQCt), row.names = colnames(dataQCt))
)

# Write to compressed h5ad file
write_h5ad(ad,paste0(datapath, "primary_MGE_SST.h5ad"),compression='gzip')

# Check file size. File MUST be <500MB to upload for MapMyCells
print(paste("Size in MB:",round(file.size(paste0(datapath, "primary_MGE_SST.h5ad"))/2^20)))

# Mapping result
mapping <- read.csv(paste0(resultdir,
                           "WangNature/MGE_SST/mapmycell/primary_MGE_SST_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1768582244791.csv"),comment.char="#")
summary(as.factor(mapping %>% filter(supercluster_bootstrapping_probability > 0.8) %>% pull (supercluster_name)))
# Amygdala excitatory 2 
# CGE interneuron 16 
# Committed oligodendrocyte precursor 7 
# Deep-layer corticothalamic and 6b 2 
# Deep-layer intratelencephalic 3 
# Deep-layer near-projecting 1 
# Eccentric medium spiny neuron 1 
# Fibroblast 2 
# LAMP5-LHX6 and Chandelier 2 
# Medium spiny neuron 2 
# MGE interneuron 8120 
# Miscellaneous 474 
# Oligodendrocyte precursor 1 
# Splatter 196 
# Upper-layer intratelencephalic 10 
mapping %>% filter(supercluster_bootstrapping_probability > 0.8) %>% dim() #8839
mapping %>% filter(cluster_bootstrapping_probability > 0.8) %>% dim() #6609
mapping_clean <- mapping %>%
  mutate(
    supercluster_name = if_else(
      supercluster_bootstrapping_probability > 0.8,
      supercluster_name,
      NA_character_
    ),
    cluster_name = if_else(
      cluster_bootstrapping_probability > 0.8,
      cluster_name,
      NA_character_
    ),
    subcluster_name = if_else(
      subcluster_bootstrapping_probability > 0.8,
      subcluster_name,
      NA_character_
    )
  )
primary_MGE_SST <- AddMetaData(primary_MGE_SST, mapping_clean)
fwrite(primary_MGE_SST@meta.data,paste0(resultdir,"WangNature/MGE_SST/mapmycell/primary_MGE_SST_meta.data.csv"), sep = "\t")

DimPlot_scCustom(subset(primary_MGE_SST, subset = !is.na(supercluster_name)), reduction = "umap", group.by = "supercluster_name", raster = FALSE, label = F)
DimPlot_scCustom(subset(primary_MGE_SST, subset = !is.na(cluster_name)), reduction = "umap", group.by = "cluster_name", raster = FALSE, label = F)
DimPlot_scCustom(subset(primary_MGE_SST, subset = subtype != "Other" & !is.na(cluster_name)), reduction = "umap", group.by = "cluster_name", raster = FALSE, label = F)
DimPlot_scCustom(subset(primary_MGE_SST, subset = (supercluster_name == "MGE interneuron" | supercluster_name == "Miscellaneous")& !is.na(cluster_name)), reduction = "umap", group.by = "cluster_name", raster = FALSE, label = F)


#pseudobulk DEG analysis ----
primary_MGE_SST_pseudo = subset(primary_MGE_SST, subset = subtype != "Other")
DefaultAssay(primary_MGE_SST_pseudo) <- "RNA"
primary_MGE_SST_pseudo[["SCT"]] <- NULL
primary_MGE_SST_pseudo[["ATAC"]] <- NULL
primary_MGE_SST_pseudo[["integrated"]] <- NULL
primary_MGE_SST_pseudo.sce <- as.SingleCellExperiment(primary_MGE_SST_pseudo)
primary_MGE_SST_pseudo.sce$ident <- NULL

#aggregate across cells by subtype and dataset(donor)
summed <- aggregateAcrossCells(primary_MGE_SST_pseudo.sce, ids = colData(primary_MGE_SST_pseudo.sce)[,c("subtype", "dataset")])

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(summed), samples=colData(summed))

#remove pseudobulk cells consisting of less than 50 cells
discarded <- summed$ncells < 50
y <- y[,!discarded]
# Another typical step in bulk RNA-seq analyses is to remove genes that are lowly expressed. This reduces computational work, improves the accuracy of mean-variance trend modelling and decreases the severity of the multiple testing correction. Here, we use the filterByExpr() function from edgeR to remove genes that are not expressed above a log-CPM threshold in a minimum number of samples (determined from the size of the smallest treatment group in the experimental design).
keep <- filterByExpr(y, group=summed$subtype)
y <- y[keep,]
summary(keep)
# Mode   FALSE    TRUE 
# logical   26397   10204

#estiamate dispersion
sizeFactors <- calcNormFactors(y$counts, method="TMM") # considerate composition bias - scaling factor like library size
disp <- setNames(edgeR::estimateDisp(y$counts)$tagwise.dispersion, rownames(y$counts)) # tagwise dispersion for every genes

#perform DE test using a generalized linear mixed model using glmmSeq package, LRT test was used here with a reduced model without the "subtype" covariate
results <- glmmSeq(modelFormula = ~ subtype + log2_age + (1 | dataset),
                   reduced = ~ log2_age + (1 | dataset),
                   countdata = y$counts,
                   metadata = y$samples,
                   dispersion = disp,
                   sizeFactors = sizeFactors,
                   cores = 12,
                   progress = TRUE)
results <- glmmQvals(results)
# LRT
# Not Significant     Significant 
# 229                  9975 
stats <- data.frame(summary(results))
colnames(stats)[c(14)] <- c("Q_LRT")
# add log2 fold change information from predicted values
stats$l2fc <- stats$subtypeIN.MGE.SST.2 / log(2) 
write.csv(stats, paste0(resultdir,
                        "WangNature/MGE_SST/glmmseq_MGE.SST.2_vs_MGE.SST.1_res.csv"))


# plot -----
## volcono plot -----
stats$p.adj <- p.adjust(stats$P_LRT, method = 'BH')

stats$significant <- ifelse(stats$p.adj < 0.05, TRUE, FALSE)
stats <- stats %>% mutate(gene_type = case_when(l2fc > 1 & p.adj < 0.05 ~ "SST-deep biased",
                                                l2fc < -1 & p.adj < 0.05 ~ "SST-upper biased",
                                                TRUE ~ "Other"))
stats$gene_type <- factor(stats$gene_type, levels = c("SST-deep biased", "SST-upper biased", "Other"))
stats$symbol <- rownames(stats)
write.csv(stats, paste0(resultdir,
                        "WangNature/MGE_SST/glmmseq_MGE.SST.2_vs_MGE.SST.1_res_final.csv"), row.names = F)

write.csv(stats %>% filter(gene_type != "Other"), 
          paste0(resultdir,"WangNature/MGE_SST/glmmseq_MGE.SST.2_vs_MGE.SST.1_res_final_sig.csv"), row.names = F)

# select genes to highlight in volcano plot
genes_to_highlight <- stats %>%
  filter((l2fc > 2.3 & p.adj < 0.005) | (l2fc < -2.3 & p.adj < 0.005))

# Add colour, size and alpha (transparency) to volcano plot 
cols <- c("SST.2 biased" = "#1d3557", "SST.1 biased" = "#e63946", "Other" = "grey")

ggplot(stats, aes(x = l2fc, y = -log10(p.adj))) +
  geom_point(mapping = aes(color = gene_type), alpha = 0.2, size = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_size_manual(values = 1) + # Modify point size
  # scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
  #                    limits = c(-5.5, 5.5)) +
  geom_point(data = genes_to_highlight,
             mapping = aes(color = gene_type),
             size = 2,
             alpha = 1) +
  scale_color_manual(values = cols) + # Modify point colour
  geom_text_repel(data = genes_to_highlight, # Add labels last to appear as the top layer  
                  aes(label = symbol),
                  size = 2,
                  nudge_x = 0.2,
                  nudge_y = 0.4,
                  max.overlaps = 20,
                  min.segment.length = 0.21
  ) +
  # ylim(0,4) +
  # xlim(-8,8) +
  labs(x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       color = "Expression \ndifference") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.text=element_text(size=15),
        #legend.title=element_text(size=18)
  )
ggsave(filename = paste0(resultdir,
                         "WangNature/MGE_SST/volcano.pdf"), width = 7, height = 5, units = "in")

## plot top markers in umap -----
DefaultAssay(primary_MGE_SST) <- "SCT"

## Inh L1-3 SST FAM20A
for (i in c( "LINGO2", "SST", "CDH12", "SYTL5", "NCAM2", "SOX6", "TAC1", "KIRREL3", "FAM20A")) {
  p <- FeaturePlot_scCustom(primary_MGE_SST, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}

## Inh L5-6 SST PAWR
for (i in c("NDST4", "LRRTM4", "NMU", "PTPRT", "DCN", "ALCAM", "KLHL14", "GALNTL6", "TMCC3", "SST","PAWR")) {
  p <- FeaturePlot_scCustom(primary_MGE_SST, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}

for (i in c("TACR1", "NPY", "ADRA1A", "PENK", "KIRREL3", "ANOS1", "AL121578.3", "SYTL5")) {
  p <- FeaturePlot_scCustom(primary_MGE_SST, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}

for (i in c("CALB1", "FRZB", "ADGRG6", "STK32A", "HPGD", 
            "B3GAT2", "NPY", "KLHDC8A", "NPM1P10", "TH", "GXYLT2", "MIR548F2", "SNCG")) {
  p <- FeaturePlot_scCustom(primary_MGE_SST, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}

# SST-1/2 marker genes top 10 
library(readr)
type_markers_SST <- fread("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/MGE_SST/type_markers.csv")
type_markers_SST$p_val_adj <- as.numeric(type_markers_SST$p_val_adj)
type_markers_SST_filt <- type_markers_SST %>% filter(avg_log2FC > 0 & p_val_adj <0.05 & pct.1 > 0.25) %>% 
  group_by(cluster) %>% arrange(-avg_log2FC) %>% slice_head(n=50)
type_markers_SST_1_filt <- type_markers_SST_filt %>% filter(cluster == "IN-MGE-SST-1") %>% arrange(avg_log2FC)
type_markers_SST_2_filt <- type_markers_SST_filt %>% filter(cluster == "IN-MGE-SST-2") %>% arrange(avg_log2FC)
features_SST1_top10 <-type_markers_SST_1_filt$gene[1:10]
features_SST2_top10 <-type_markers_SST_2_filt$gene[1:10]

for (i in features_SST1_top10) {
  p <- FeaturePlot_scCustom(primary_MGE_SST, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}

for (i in features_SST2_top10) {
  p <- FeaturePlot_scCustom(primary_MGE_SST, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}

# combine plot
primary_MGE_SST_ageGroup_plot <- DimPlot_scCustom(primary_MGE_SST, reduction = "umap", group.by = "Group", raster = F,colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8")) + NoLegend() + NoAxes() + labs(title = NULL)
primary_MGE_SST_subtype_plot <- DimPlot_scCustom(primary_MGE_SST, reduction = "umap", group.by = "subtype", raster = F, colors_use = c("#e63946", "#1d3557","grey")) + NoLegend() + NoAxes() + labs(title = NULL)
primary_MGE_SST_region_plot <- DimPlot_scCustom(primary_MGE_SST, reduction = "umap", group.by = "region_summary", raster = F, colors_use = c("#009e73", "#ffa500", "#0072b2")) + NoLegend() + NoAxes() + labs(title = NULL)
primary_MGE_SST_RVAS_scDRS_plot <- FeaturePlot_scCustom(primary_MGE_SST, reduction = "umap", 
                                                        features = paste0("scDRS_run_v8_253_", trait, "_norm_score"),colors_use = c("#053061", "#2166AC", "#4393C3", "#92C5DE","#F7F7F7",
                                                                                                                                                             "#F4A582", "#D6604D", "#B2182B", "#67001F"),max.cutoff = 5,min.cutoff = -3,na_cutoff  = NULL,raster = F) + NoLegend() + NoAxes() + labs(title = NULL)
primary_MGE_SST_GWAS_scDRS_plot <- FeaturePlot_scCustom(primary_MGE_SST, reduction = "umap", 
                                                        features = "scDRS_run_GWAS_v3_100_ASD_Matoba_2020_norm_score",colors_use = c("#053061", "#2166AC", "#4393C3", "#92C5DE","#F7F7F7",
                                                                                                                                                     "#F4A582", "#D6604D", "#B2182B", "#67001F"),max.cutoff = 5,min.cutoff = -3,na_cutoff  = NULL,raster = F) + NoLegend() + NoAxes() + labs(title = NULL)

primary_MGE_SST_ageGroup_plot + primary_MGE_SST_region_plot + primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + primary_MGE_SST_GWAS_scDRS_plot + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_combined_",trait,".pdf"), width = 7.5, height = 5, units = "in")

primary_MGE_SST_ageGroup_plot + primary_MGE_SST_region_plot + primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + TACR1 + NPY + ADRA1A + PENK + KIRREL3 + ANOS1 + AL121578.3 + SYTL5 + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_markers_combined_",trait,".pdf"), width = 7.5, height = 10, units = "in")

primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + 
  NPY + FRZB + B3GAT2 + STK32A + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_markers_combined_selected_",trait,".pdf"), width = 7.5, height = 5, units = "in")

primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + 
  SST + FAM20A + KIRREL3 + SYTL5 + LINGO2 + CDH12  + NCAM2 + SOX6 + TAC1 + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_markers_combined_",trait,"_L1-3-SST-FAM20A.pdf"), width = 7.5, height = 10, units = "in")

primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + 
  SST + PAWR + ALCAM + GALNTL6 + NDST4 + LRRTM4  + NMU + PTPRT + DCN + KLHL14 + TMCC3 + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_markers_combined_",trait,"_L5-6-SST-PAWR.pdf"), width = 7.5, height = 12.5, units = "in")

primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + CDH12 + TACR1 + IGF1 + CDH7 + RBMS3 + MAN1A1 + DKK2 + ZNF385D + VWC2 + AC079296.1 + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_markers_combined_",trait,"_SST1_top10marker.pdf"), width = 7.5, height = 10, units = "in")

primary_MGE_SST_subtype_plot + primary_MGE_SST_RVAS_scDRS_plot + WLS + `RBMS3-AS3` + TAC1 + STXBP6 + PRLR + GALNTL6 + `TRHDE-AS1` + PLCH1 + GAD1 + ARX + plot_layout(ncol = 3)
ggsave(file=paste0(resultdir,
                   "WangNature/MGE_SST/subtype_markers_combined_",trait,"_SST2_top10marker.pdf"), width = 7.5, height = 10, units = "in")


