# ── Load data ──────────────────────────────────────────────────────────────
count_ddid   <- read_excel("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/counts_asd_ddid.xlsx")
count_noddid <- read_excel("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/counts_asd_noddid.xlsx")
kaplanis <- fread("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/Kaplanis_DDD_ASD_counts_by_gene_new_mis_cats_2025-04-09.txt") %>%
  transmute(Gene,
            PTV_Proband  = PTV_Proband_Male  + PTV_Proband_Female,
            Mis2_Proband = Mis2_Proband_Male + Mis2_Proband_Female,
            Mis1_Proband = Mis1_Proband_Male + Mis1_Proband_Female)


rvas_fdr001 <- fread(file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001.txt")
rvas_fdr01 <- fread(file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr01.txt")
RVAS_result = readRDS("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/2026.01.07_ASD-Zscore-FDRtransform.RDS")

count_ddid <- count_ddid %>%
  left_join(kaplanis, by = "Gene", suffix = c("", "_kap")) %>%
  mutate(
    PTV_Proband  = PTV_Proband  + coalesce(PTV_Proband_kap,  0L),
    Mis2_Proband = Mis2_Proband + coalesce(Mis2_Proband_kap, 0L),
    Mis1_Proband = Mis1_Proband + coalesce(Mis1_Proband_kap, 0L)
  ) %>%
  select(-ends_with("_kap"))

# Filter to FDR < 0.01 genes
count_ddid_fdr01   <- count_ddid   %>% filter(Gene %in% rvas_fdr01$gene)
count_noddid_fdr01 <- count_noddid %>% filter(Gene %in% rvas_fdr01$gene)

# Filter to FDR < 0.001 genes
count_ddid   <- count_ddid   %>% filter(Gene %in% rvas_fdr001$gene)
count_noddid <- count_noddid %>% filter(Gene %in% rvas_fdr001$gene)

# ── Analysis: Per-proband mutation rate ──────────────────────────────────────
# compare per-proband mutation rates between DDID and NoDDID cohorts
library(purrr)
# Cohort sizes
n_ddid_probands   <- 13841
n_noddid_probands <- 24839

classify_perproband <- function(count_ddid, count_noddid, cols, col_suffix) {
  col_ddid   <- paste0("count_", col_suffix, "_ddid")
  col_noddid <- paste0("count_", col_suffix, "_noddid")
  
  df <- merge(
    count_ddid   %>% mutate(!!col_ddid   := rowSums(across(all_of(cols)))) %>% select(Gene, all_of(col_ddid)),
    count_noddid %>% mutate(!!col_noddid := rowSums(across(all_of(cols)))) %>% select(Gene, all_of(col_noddid))
  ) %>%
    filter(.data[[col_ddid]] + .data[[col_noddid]] > 0) %>%  
    mutate(
      # Per-proband mutation rates
      rate_ddid   = .data[[col_ddid]]   / n_ddid_probands,
      rate_noddid = .data[[col_noddid]] / n_noddid_probands,
      # p_hat for reference
      p_hat       = .data[[col_ddid]] / (.data[[col_ddid]] + .data[[col_noddid]]),
      # Classification based on per-proband rate
      Group_perproband = ifelse(rate_ddid > rate_noddid, "DDID", "NoDDID")
    )
  
  message("Per-proband group counts (", col_suffix, "):")
  print(table(df$Group_perproband))
  
  df
}

# ── Save gene lists ────────────────────────────
save_ddid_gene_lists_custom <- function(df, group_col, RVAS_result, datadir, suffix) {
  for (grp in c("DDID", "NoDDID")) {
    genes <- df$Gene[df[[group_col]] == grp]
    out <- RVAS_result %>%
      filter(gene_name %in% genes) %>%
      select(gene_name, p_fdr) %>%
      arrange(p_fdr) %>%
      rename(GENE = gene_name,
             !!paste0("ASC_Pfdr_", tolower(grp), "_", suffix) := p_fdr)
    fwrite(
      out,
      file = paste0(datadir, "RVAS_result/full_results_wcounts_pfdr_",
                    tolower(grp), "_", suffix, ".tsv"),
      sep  = "\t"
    )
  }
  message("Saved gene lists for: ", suffix)
}

count_perproband <- classify_perproband(
  count_ddid, count_noddid,
  cols       = c("PTV_Proband", "Mis2_Proband", "Del_proband"),
  col_suffix = "new_DMN"
)
# Per-proband group counts (new_DMN):
#   
# DDID NoDDID 
# 205     44 
fwrite(count_perproband,file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_perproband_group.txt",sep  = "\t")

count_perproband_fdr01 <- classify_perproband(
  count_ddid_fdr01, count_noddid_fdr01,
  cols       = c("PTV_Proband", "Mis2_Proband", "Del_proband"),
  col_suffix = "new_DMN_fdr01"
)
# Per-proband group counts (new_DMN):
#   
# DDID NoDDID  
# 293    105 
fwrite(count_perproband_fdr01,file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr01_PtvMis2Del_perproband_group.txt",sep  = "\t")

save_ddid_gene_lists_custom(
  count_perproband, "Group_perproband", RVAS_result, datadir,
  suffix = "perproband_new_DMN"
)
save_ddid_gene_lists_custom(
  count_perproband_fdr01, "Group_perproband", RVAS_result, datadir,
  suffix = "perproband_new_DMN_fdr01"
)

# ── Analysis: p_hat ──────────────────────────────────────────────
# ── Helper: group genes into DDID / NoDDID based on variant count ratio ────────
#
# @param count_ddid    data.frame with ddid counts
# @param count_noddid  data.frame with noddid counts
# @param cols          character vector of count columns to sum (e.g. PTV, Mis1, Mis2)
# @param col_suffix    suffix for the new count columns (e.g. "DMN", "new_DMN")
#
# @return data.frame with Gene, count columns, p_hat, Group
classify_ddid_group <- function(count_ddid, count_noddid, cols, col_suffix) {
  col_ddid   <- paste0("count_", col_suffix, "_ddid")
  col_noddid <- paste0("count_", col_suffix, "_noddid")
  df <- merge(count_ddid %>% mutate(!!col_ddid   := rowSums(across(all_of(cols)))) %>% select(Gene, all_of(col_ddid)),
              count_noddid %>% mutate(!!col_noddid := rowSums(across(all_of(cols)))) %>% select(Gene, all_of(col_noddid))
  ) %>% mutate(p_hat = .data[[col_ddid]] / (.data[[col_ddid]] + .data[[col_noddid]]))
  p0 <- sum(df[[col_ddid]]) / (sum(df[[col_ddid]]) + sum(df[[col_noddid]]))
  message("p0 (", col_suffix, ") = ", round(p0, 4))
  df %>% mutate(Group = ifelse(p_hat > p0, "DDID", "NoDDID"))
}

count_DMN <- classify_ddid_group( count_ddid, count_noddid,
                                  cols       = c("PTV_Proband", "Mis2_Proband", "Del_proband"),
                                  col_suffix = "new_DMN")
# p0 (new_DMN) = 0.6437
# Group counts (new_DMN):
# 
# DDID NoDDID 
# 121    128  

fwrite(count_DMN,file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_group.txt",sep  = "\t")

save_ddid_gene_lists_custom(
  count_DMN, "Group", RVAS_result, datadir,
  suffix = "new_DMN"
)

