library(dplyr)
library(data.table)
library(readxl)

datadir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/"
rvas_result <- fread(paste0(datadir, "RVAS_result/full_results_wcounts_2025-10-08.txt"))
RVAS_result = readRDS("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/2026.01.07_ASD-Zscore-FDRtransform.RDS")

rvas_fdr001 = rvas_result %>% filter(FDR < 0.001 & !FLAG) %>% select(gene)
rvas_fdr01 = rvas_result %>% filter(FDR < 0.01 & !FLAG) %>% select(gene)

rvas_fdr001_flag = rvas_result %>% filter(FDR < 0.001 & FLAG) %>% select(gene)

rvas_result_pfdr <- RVAS_result %>%
  select(gene_name, p_fdr) %>%
  arrange(p_fdr) %>%
  rename(GENE = gene_name,
         ASC_Pfdr = p_fdr) %>%
  filter(!GENE %in% rvas_fdr001_flag$gene) 

fwrite(rvas_result_pfdr, file = paste0(datadir, "RVAS_result/full_results_wcounts_pfdr.tsv"), sep = "\t")
fwrite(rvas_fdr001, file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001.txt", sep = "\t")
fwrite(rvas_fdr01, file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr01.txt", sep = "\t")

