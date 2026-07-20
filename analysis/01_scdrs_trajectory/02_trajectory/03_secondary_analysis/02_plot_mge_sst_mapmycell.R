library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(readxl)
library(ggrepel)
library(tidyverse)

#import  data
datapath = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/"
resultdir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
trait = "ASC_Pfdr" 
primary_MGE_SST <- readRDS(paste0(datapath, "primary_MGE_SST.rds"))
meta <- fread(paste0(resultdir,"WangNature/MGE_SST/mapmycell/primary_MGE_SST_meta.data.csv"), na.strings = "")
science_add7046_table_s3 <- read_excel("/mnt/isilon/gandal_lab/liaoyd/data/allen_brain_atlas/WHB_taxonomy/science.add7046_table_s3.xlsx")

meta2 <- meta %>% filter((subtype == "IN-MGE-SST-1" | subtype == "IN-MGE-SST-2") & !is.na(cluster_name))
cluster_counts <- meta2 %>%
  group_by(subtype, cluster_name) %>%
  tally() %>%
  group_by(subtype) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

main_clusters <- meta2 %>%
  dplyr::count(cluster_name) %>%
  mutate(total_prop = n / sum(n)) %>%
  filter(total_prop >= 0.02) %>%
  pull(cluster_name)

plot_df <- cluster_counts %>%
  mutate(cluster_plot = ifelse(cluster_name %in% main_clusters, cluster_name, "Others")) %>%
  filter(cluster_plot != "Others") %>%
  group_by(subtype, cluster_plot) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(subtype) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

enrichment_order <- plot_df %>%
  filter(subtype %in% c("IN-MGE-SST-1", "IN-MGE-SST-2")) %>% 
  select(subtype, cluster_plot, prop) %>%
  tidyr::pivot_wider(names_from = subtype, values_from = prop, values_fill = 0) %>%
  mutate(log2_fc = log2((`IN-MGE-SST-1` + 0.001) / (`IN-MGE-SST-2` + 0.001))) %>%
  arrange(log2_fc) 

plot_df$cluster_plot <- factor(plot_df$cluster_plot, levels = enrichment_order$cluster_plot)

freq_plot <- ggplot(plot_df, aes(x = cluster_plot, y = prop, fill = subtype)) +
  geom_col(position = "fill", width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(
    values = c("IN-MGE-SST-1" = "#e63946", "IN-MGE-SST-2" = "#1d3557")
  ) +
  labs(
    title = "Cluster Specificity by SST Subtype",
    x = "10x Whole Human Brain Clusters",
    y = "Percentage Composition",
    fill = "Subtype"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

print(freq_plot)


# top 5 clusters
group1_clusters <- enrichment_order %>% arrange(desc(log2_fc)) %>% slice_head(n=5) %>% pull(cluster_plot)
group2_clusters <- enrichment_order %>% arrange(log2_fc) %>% slice_head(n=5) %>% pull(cluster_plot)

score_colors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", 
                  "#F4A582", "#D6604D", "#B2182B", "#67001F")
score_limits <- c(-3, 5)
score_feature <- paste0("scDRS_run_v8_253_", trait, "_norm_score")

unified_scale <- scale_color_gradientn(
  colors = score_colors,
  limits = score_limits,
  oob = scales::squish  
)

# Top 5 enriched clusters
p1a <- DimPlot_scCustom(
  subset(primary_MGE_SST, subset = cluster_name %in% group1_clusters),
  reduction = "umap",
  colors_use = c("#3B00FBFF", "#FE00FAFF", "#1CBE4FFF", "#3283FEFF", "#FEAF16FF"),
  group.by = "cluster_name", raster = FALSE, label = FALSE
) + NoAxes() + labs(title = "10x Human Brain Mapping cluster")

p2a <- FeaturePlot_scCustom(
  subset(primary_MGE_SST, subset = cluster_name %in% group1_clusters),
  reduction = "umap",
  features = score_feature,
  colors_use = score_colors,
  max.cutoff = 5, min.cutoff = -3, na_cutoff = NULL, raster = FALSE
) + NoAxes() + labs(title = "scDRS Mapping") +
  unified_scale  

# Top 5 depleted clusters
p1b <- DimPlot_scCustom(
  subset(primary_MGE_SST, subset = cluster_name %in% group2_clusters),
  reduction = "umap",
  colors_use = c("#B10DA1FF", "#FBE426FF", "#683B79FF", "#90AD1CFF", "#2ED9FFFF"),
  group.by = "cluster_name", raster = FALSE, label = FALSE
) + NoAxes() + labs(title = "10x Human Brain Mapping cluster")

p2b <- FeaturePlot_scCustom(
  subset(primary_MGE_SST, subset = cluster_name %in% group2_clusters),
  reduction = "umap",
  features = score_feature,
  colors_use = score_colors,
  max.cutoff = 5, min.cutoff = -3, na_cutoff = NULL, raster = FALSE
) + NoAxes() + labs(title = "scDRS Mapping") +
  unified_scale  

p1a + p2a
p1b + p2b


# annotation
feature_summary <- science_add7046_table_s3 %>%
  mutate(group = case_when(
    `Cluster name` %in% group1_clusters ~ "SST-1 (High)",
    `Cluster name` %in% group2_clusters ~ "SST-2 (Low)",
    TRUE ~ "Others"
  )) %>%
  filter(group != "Others")

fwrite(feature_summary,paste0(resultdir,"WangNature/MGE_SST/mapmycell/IN_MGE_SST_mapping.xlsx"), sep = "\t")


get_unique_peps <- function(pep_strings) {
  pep_strings %>%
    str_split(" ") %>%       
    unlist() %>%             
    unique() %>%             
    setdiff(c("", NA))       
}
get_gene_freq <- function(gene_strings) {
  if(length(gene_strings) == 0) return(data.frame(Gene=character(), Count=numeric()))
  gene_strings %>%
    str_split(",") %>% 
    unlist() %>%
    str_trim() %>%           
    table() %>%              
    as.data.frame() %>%      
    setNames(c("Gene", "Count")) %>% 
    arrange(desc(Count))
}


sst1_pep_all <- feature_summary %>% 
  filter(group == "SST-1 (High)") %>% 
  pull(`Neuropeptide auto-annotation`) %>% 
  get_unique_peps()

sst2_pep_all <- feature_summary %>% 
  filter(group == "SST-2 (Low)") %>% 
  pull(`Neuropeptide auto-annotation`) %>% 
  get_unique_peps()


common_peptides <- intersect(sst1_pep_all, sst2_pep_all)
sst1_only_peptides <- setdiff(sst1_pep_all, sst2_pep_all)
sst2_only_peptides <- setdiff(sst2_pep_all, sst1_pep_all)

sst1_gene_table <- feature_summary %>% 
  filter(group == "SST-1 (High)") %>% 
  pull(`Top Enriched Genes`) %>% 
  get_gene_freq()

sst2_gene_table <- feature_summary %>% 
  filter(group == "SST-2 (Low)") %>% 
  pull(`Top Enriched Genes`) %>% 
  get_gene_freq()

sst1_distinctive_genes <- sst1_gene_table %>%
  filter(Count >= 2) %>%
  filter(!Gene %in% sst2_gene_table$Gene[sst2_gene_table$Count >= 1])

sst2_distinctive_genes <- sst2_gene_table %>%
  filter(Count >= 2) %>%
  filter(!Gene %in% sst1_gene_table$Gene[sst1_gene_table$Count >= 1])

common_high_freq <- inner_join(
  sst1_gene_table %>% filter(Count >= 2),
  sst2_gene_table %>% filter(Count >= 2),
  by = "Gene", suffix = c("_SST1", "_SST2")
)

sst1_only_peptides
sst2_only_peptides
common_peptides
common_high_freq
sst1_distinctive_genes
# Inh L1-2 SST CCNJL
# CMTM8/MYO5B

#MEF2C；PENK;CMTM8/MYO5B
FeaturePlot_scCustom(primary_MGE_SST, features = "PENK", reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/PENK", ".pdf"), width = 5, height = 5, units = "in")

FeaturePlot_scCustom(primary_MGE_SST, features = "NPPC", reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/NPPC", ".pdf"), width = 5, height = 5, units = "in")

FeaturePlot_scCustom(primary_MGE_SST, features = "MEF2C", reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/MEF2C", ".pdf"), width = 5, height = 5, units = "in")

FeaturePlot_scCustom(primary_MGE_SST, features = "MYO5B", reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/MYO5B", ".pdf"), width = 5, height = 5, units = "in")

FeaturePlot_scCustom(primary_MGE_SST, features = "SLC17A8", reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = paste0(resultdir,"WangNature/MGE_SST/markers/SLC17A8", ".pdf"), width = 5, height = 5, units = "in")

