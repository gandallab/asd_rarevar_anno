library(Seurat)
library(scCustomize)
library(ggplot2)
library(data.table)
library(dplyr)

datapath = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"

# A Seurat object containing count matrices, dimension reduction, and metadata  from 38 samples in this experiment (snMultiome_atlas_Seurat_object.rds) was downloaded from https://datadryad.org/dataset/doi:10.5061/dryad.2280gb612#citations.
obj = readRDS(paste0(datapath, "snMultiome_atlas_Seurat_object.rds")) # 74GB

gene_numbers <- c(253,416,696,951)
traits = c("ASC_Pfdr")
version = "run_v8"
for (gene_number in gene_numbers) {
  for (trait in traits){
score_path <- paste0(
  "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/scDRS/",
  version, "/score_file/",
  gene_number, "_genes/", trait, ".full_score.gz"
)
score_dt <- fread(score_path)
score_df <- as.data.frame(score_dt)
rownames(score_df) <- score_df$V1

scDRS_norm_score <- score_df %>%
  dplyr::select(norm_score, nlog10_pval)
score_col <- paste0("scDRS_", version, "_", gene_number, "_", trait, "_norm_score")
pval_col  <- paste0("scDRS_", version, "_", gene_number, "_", trait, "_nlog10_pval")

obj <- AddMetaData(
  obj,
  metadata = scDRS_norm_score,
  col.name = c(score_col, pval_col)
)

p_scDRS <- FeaturePlot_scCustom(
  obj,
  features   = score_col,
  reduction  = "wnn.umap",
  colors_use = c("#2166AC", "#4393C3", "#92C5DE", "#CADDED",
                 "#FDDBC7",
                 "#F4A582", "#D6604D", "#B2182B", "#67001F"),
  max.cutoff = 5,
  min.cutoff = -3,
  na_cutoff  = NULL,
  raster     = TRUE 
)
p_scDRS_build <- ggplot_build(p_scDRS)
df_umap <- p_scDRS_build$data[[1]]

pdf(
  file = paste0(resultdir,
                "WangNature/scDRS/",version,"/",version,"_norm_score_umap_",
                gene_number, "genes_",trait,".pdf"),
  width  = 2.4,
  height = 2.9
  # units  = "in",
  # res    = 300
)
plot(df_umap$x, df_umap$y,
     xlab = "",
     ylab = "",
     col  = scales::alpha(df_umap$colour,0.5),
     asp  = NA,
     pch  = 16,
     cex  = 0.05,
     axes = FALSE,
     xaxs = "i",
     yaxs = "i")
dev.off()
  }
}

# subclass umap
DimPlot_scCustom(obj, group.by = "subclass", reduction = "wnn.umap", pt.size = 1, raster = T)
ggsave(file = paste0(resultdir,"WangNature/subclass_umap_filtered.pdf"), width = 8, height = 6, units = "in")

# plot the scDRS associated Z-score ----
traits = c("ASC_Pfdr")
for (trait in traits){
scDRS_RVAS_sex_fdr001 <- fread(file = paste0(
  resultdir,"WangNature/scDRS/run_v8/downstream_analysis/253_genes/", trait, ".scdrs_group.subclass")) %>%
  mutate(assoc_mcp_fdr = p.adjust(assoc_mcp,method = "fdr"),
         prop_fdr_05 = 100*n_fdr_0.05/n_cell,
         prop_fdr_1 = 100*n_fdr_0.1/n_cell,
         prop_fdr_2 = 100*n_fdr_0.2/n_cell)
fwrite(scDRS_RVAS_sex_fdr001, file = paste0(resultdir,"WangNature/scDRS/run_v8/downstream_analysis/253_genes/", trait, ".scdrs_group.subclass.csv"), sep = "\t")

df= scDRS_RVAS_sex_fdr001 %>% select(group,assoc_mcz,assoc_mcp_fdr) %>% arrange(-assoc_mcz) %>%
  mutate(group = factor(group, levels = group),
                        sig = case_when(
                          assoc_mcp_fdr < 0.001 ~ "***",
                          assoc_mcp_fdr < 0.01  ~ "**",
                          assoc_mcp_fdr < 0.05  ~ "*",
                          TRUE ~ "")
                        )

z_lim = range(scDRS_RVAS_sex_fdr001$assoc_mcz)
ggplot(df, aes(x = group, y = "assoc_mcz", fill = assoc_mcz)) +
  geom_tile(width = 0.9, height = 0.9) +
  geom_text(aes(label = sig),
    color = "white",
    size = 3) +
  scale_fill_gradientn(colors = 
    c("#08306B",     
      "#6BAED6", 
      "#C6DBEF",   
      "white",
      "#F9C6A3",
      "#F07020",  
      "#E03060", 
      "#C0267A"),
    
    values = scales::rescale(c(z_lim[1], -1, -0.5, 0, 1.3, 2, 3.5, z_lim[2])),
    limits = z_lim
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = paste0(resultdir,
                         "WangNature/scDRS/run_v8/scDRS_RVAS_sex_fdr001_zscore_",trait,".pdf"), width = 7, height = 2)

}
