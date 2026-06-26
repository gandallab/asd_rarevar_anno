library(data.table)
library(dplyr)
library(tidyr)

# ── Column map: only _sex score types ─────────────────────────────────────────
make_colmap <- function(trait) {
  tibble::tibble(
    score_type = c(
      "GWAS_sex",
      "RVAS_sex_fdr001",
      "RVAS_sex_fdr01",
      "RVAS_sex_fdr05",
      "RVAS_sex_fdr1"),
    colname = c(
      "scDRS_run_GWAS_v3_norm_score",
      paste0("scDRS_run_v8_253_",              trait, "_norm_score"),
      paste0("scDRS_run_v8_416_",              trait, "_norm_score"),
      paste0("scDRS_run_v8_696_",              trait, "_norm_score"),
      paste0("scDRS_run_v8_951_",              trait, "_norm_score")
    )
  )
}


# ── Helper: process one cell type (EN or IN) ──────────────────────────────────
#
# @param cell_type    "EN" or "IN"
# @param col_score1   column name in metadata for score1 (e.g. GWAS)
# @param col_score2   column name in metadata for score2 (e.g. RVAS)
# @param score1_label score_type label for score1 (e.g. "GWAS_sex")
# @param score2_label score_type label for score2 (e.g. "RVAS_sex_fdr001")
# @param version1     scDRS run version for score1 bin50 file path
# @param version2     scDRS run version for score2 bin50 file path
# @param gene_number1 gene set size or label for score1 bin50 file path
# @param gene_number2 gene set size or label for score2 bin50 file path
# @param trait1       trait identifier for score1 bin50 file name
# @param trait2       trait identifier for score2 bin50 file name
# @param resultdir    base path to result directory
#
# @return named list: bin50_scDRS, scDRS_long, tt_res,
#                     neglog10_fdr_enrich, neglog10_fdr_comp, df_v1_age

process_cell_type <- function(cell_type,
                              col_score1,    col_score2,
                              score1_label,  score2_label,
                              version1,      version2,
                              gene_number1,  gene_number2,
                              trait1,        trait2,
                              resultdir) {
  
  # ── 1. Bin-level enrichment FDR ─────────────────────────────────────────────
  
  read_bin50 <- function(version, gene_number, trait, label) {
    fread(file.path(
      resultdir,
      "WangNature/scDRS", version, "downstream_analysis",
      paste0(gene_number, "_genes"), cell_type,
      paste0(trait, ".scdrs_group.pseudotime_bin50")
    )) %>%
      mutate(
        assoc_mcp_fdr = p.adjust(assoc_mcp, method = "fdr"),
        score_type    = label
      )
  }
  
  bin50_score1 <- read_bin50(version1, gene_number1, trait1, score1_label)
  bin50_score2 <- read_bin50(version2, gene_number2, trait2, score2_label)
  
  bin50_scDRS <- rbind(bin50_score1, bin50_score2) %>%
    group_by(score_type) %>%
    mutate(bin_id = group + 1) %>%
    ungroup()
  
  neglog10_fdr_enrich <- -log10(pmax(bin50_scDRS$assoc_mcp_fdr, 1e-300))
  
  # ── 2. Load metadata and pseudotime bin edges ────────────────────────────────
  
  metadata <- fread(file.path(
    resultdir, "WangNature/slingshot", cell_type,
    paste0(cell_type, "_lineage_slingshot_results.csv")
  )) %>% as.data.frame()
  
  message("col_score1 found: ", col_score1 %in% names(metadata))
  message("col_score2 found: ", col_score2 %in% names(metadata))
  
  bin_edges <- fread(file.path(
    resultdir, "WangNature/slingshot", cell_type,
    paste0(cell_type, "_pseudotime_edges_50.txt")
  ))
  
  breaks_vec     <- bin_edges$V1
  breaks_vec[1]  <- breaks_vec[1]  - 1e-6
  breaks_vec[51] <- breaks_vec[51] + 1e-6
  
  # ── 3. Reshape to long format ────────────────────────────────────────────────
  
  col_map_pair <- tibble::tibble(
    score_type = c(score1_label, score2_label),
    colname    = c(col_score1,   col_score2)
  ) %>%
    filter(colname %in% names(metadata))
  
  rename_vec <- setNames(col_map_pair$score_type, col_map_pair$colname)
  
  scDRS_long <- metadata %>%
    select(ID, average_pseudotime, log2_age, all_of(col_map_pair$colname)) %>%
    rename_with(.fn = ~ rename_vec[.x], .cols = all_of(names(rename_vec))) %>%
    pivot_longer(
      cols      = -c(ID, average_pseudotime, log2_age),
      names_to  = "score_type",
      values_to = "score"
    ) %>%
    mutate(
      bin_id = cut(average_pseudotime,
                   breaks = breaks_vec,
                   include.lowest = TRUE,
                   labels = FALSE)
    )
  
  # ── 4. Pseudobulk paired t-test per bin ─────────────────────────────────────
  
  scDRS_pb <- metadata %>%
    mutate(
      bin_id = cut(average_pseudotime,
                   breaks = breaks_vec,
                   include.lowest = TRUE,
                   labels = FALSE)
    ) %>%
    group_by(dataset, bin_id) %>%
    summarise(
      score1   = mean(.data[[col_score1]], na.rm = TRUE),
      score2   = mean(.data[[col_score2]], na.rm = TRUE),
      log2_age = mean(log2_age,            na.rm = TRUE),
      n_cells  = n(),
      .groups  = "drop"
    ) %>%
    filter(n_cells >= 10)
  
  tt_res <- scDRS_pb %>%
    group_by(bin_id) %>%
    summarise(
      n_dataset      = n(),
      mean_diff      = mean(score1 - score2, na.rm = TRUE),
      ttest          = list(t.test(score1, score2, paired = TRUE)),
      .groups        = "drop"
    ) %>%
    mutate(
      t_stat           = sapply(ttest, function(x) x$statistic),
      p_bin_compared   = sapply(ttest, function(x) x$p.value),
      fdr_bin_compared = p.adjust(p_bin_compared, method = "fdr")
    )
  
  neglog10_fdr_comp <- -log10(pmax(tt_res$fdr_bin_compared, 1e-300))
  
  # ── 5. Per-bin V1 proportion and mean developmental age ─────────────────────
  
  df_v1_age <- metadata %>%
    filter(!is.na(average_pseudotime)) %>%
    mutate(
      bin_id = cut(average_pseudotime,
                   breaks = breaks_vec,
                   include.lowest = TRUE,
                   labels = FALSE)
    ) %>%
    group_by(bin_id) %>%
    summarise(
      prop_V1  = mean(region_summary == "V1", na.rm = TRUE),
      age_mean = mean(log2_age, na.rm = TRUE),
      .groups  = "drop"
    )
  
  list(
    bin50_scDRS         = bin50_scDRS,
    scDRS_long          = scDRS_long,
    tt_res              = tt_res,
    neglog10_fdr_enrich = neglog10_fdr_enrich,
    neglog10_fdr_comp   = neglog10_fdr_comp,
    df_v1_age           = df_v1_age
  )
}



# ── Main function: run both cell types for one comparison and save .RData ──────
#
# @param cfg        one entry from the comparisons config list (see below)
# @param resultdir  base path to result directory
run_comparison <- function(cfg, resultdir) {
  
  message("Running: ", cfg$out_label)
  
  args <- list(
    col_score1   = cfg$col_score1,  col_score2   = cfg$col_score2,
    score1_label = cfg$score1_label, score2_label = cfg$score2_label,
    version1     = cfg$version1,    version2     = cfg$version2,
    gene_number1 = cfg$gene_number1, gene_number2 = cfg$gene_number2,
    trait1       = cfg$trait1,      trait2       = cfg$trait2,
    resultdir    = resultdir
  )
  
  EN <- do.call(process_cell_type, c(list(cell_type = "EN"), args))
  IN <- do.call(process_cell_type, c(list(cell_type = "IN"), args))
  
  range_fdr_global <- range(
    c(EN$neglog10_fdr_enrich, IN$neglog10_fdr_enrich), na.rm = TRUE)
  range_fdr_comp_global <- range(
    c(EN$neglog10_fdr_comp,   IN$neglog10_fdr_comp),   na.rm = TRUE)
  range_prop_V1_global <- range(
    c(EN$df_v1_age$prop_V1,   IN$df_v1_age$prop_V1),   na.rm = TRUE)
  range_age_mean_global <- range(
    c(EN$df_v1_age$age_mean,  IN$df_v1_age$age_mean),  na.rm = TRUE)
  
  scDRS_EN_long  <- EN$scDRS_long;  scDRS_IN_long  <- IN$scDRS_long
  bin50_EN_scDRS <- EN$bin50_scDRS; bin50_IN_scDRS <- IN$bin50_scDRS
  tt_res_EN      <- EN$tt_res;      tt_res_IN      <- IN$tt_res
  df_v1_age_EN   <- EN$df_v1_age;   df_v1_age_IN   <- IN$df_v1_age
  
  out_rdata <- file.path(
    resultdir, "WangNature/slingshot",
    paste0("scDRS_RVAS_GWAS_norm_score_", cfg$out_label, ".RData")
  )
  
  save(range_fdr_global, range_fdr_comp_global,
       range_prop_V1_global, range_age_mean_global,
       scDRS_EN_long, scDRS_IN_long,
       bin50_EN_scDRS, bin50_IN_scDRS,
       tt_res_EN, tt_res_IN,
       df_v1_age_EN, df_v1_age_IN,
       file = out_rdata)
  message("Saved: ", out_rdata)
  
  fwrite(bin50_EN_scDRS,
         file = file.path(resultdir, "WangNature/slingshot",
                          paste0("scDRS_RVAS_GWAS_norm_score_EN_", cfg$out_label, ".csv")),
         sep = "\t")
  fwrite(bin50_IN_scDRS,
         file = file.path(resultdir, "WangNature/slingshot",
                          paste0("scDRS_RVAS_GWAS_norm_score_IN_", cfg$out_label, ".csv")),
         sep = "\t")
  
  invisible(list(EN = EN, IN = IN,
                 range_fdr_global      = range_fdr_global,
                 range_fdr_comp_global = range_fdr_comp_global,
                 range_prop_V1_global  = range_prop_V1_global,
                 range_age_mean_global = range_age_mean_global))
}



# ── Comparison config list ─────────────────────────────────────────────────────
#
# Each entry defines one pairwise comparison. To add a new comparison,
# simply append a new list() entry here — no function changes needed.
#
# Fields:
#   col_score1/2    : exact column names in the metadata CSV
#   score1/2_label  : score_type values that appear in scDRS_long and bin50_scDRS
#                     (must match what your plotting code uses in scale_color_manual)
#   version1/2      : scDRS run version for each score (for bin50 path)
#   gene_number1/2  : gene set size or label for each score (for bin50 path)
#   trait1/2        : trait identifier for each score bin50 filename
#   out_label       : used for .RData and .csv filenames
#
# bin50 path template (applied independently to score1 and score2):
#   {resultdir}/WangNature/scDRS/{version}/downstream_analysis/
#   {gene_number}_genes/{cell_type}/{trait}.scdrs_group.pseudotime_bin50

GWAS_Matoba_100  <- list(version = "run_GWAS_v3", gene_number = 100,  trait = "ASD_Matoba_2020")
GWAS_Matoba_200  <- list(version = "run_GWAS_v3", gene_number = 200,  trait = "ASD_Matoba_2020")
RVAS_fdr001      <- list(version = "run_v8",      gene_number = 253,  trait = "ASC_Pfdr")
RVAS_fdr01       <- list(version = "run_v8",      gene_number = 416,  trait = "ASC_Pfdr")
RVAS_fdr05       <- list(version = "run_v8",      gene_number = 696,  trait = "ASC_Pfdr")
RVAS_fdr1        <- list(version = "run_v8",      gene_number = 951,  trait = "ASC_Pfdr")

make_cfg <- function(s1, s2, col1, col2, label1, label2, out_label) {
  list(
    col_score1   = col1,        col_score2   = col2,
    score1_label = label1,      score2_label = label2,
    version1     = s1$version,  version2     = s2$version,
    gene_number1 = s1$gene_number, gene_number2 = s2$gene_number,
    trait1       = s1$trait,    trait2       = s2$trait,
    out_label    = out_label
  )
}

comparisons <- list(
  
  # GWAS (Matoba 100) vs RVAS FDR thresholds
  make_cfg(GWAS_Matoba_100, RVAS_fdr001,
           "scDRS_run_GWAS_v3_100_ASD_Matoba_2020_norm_score",
           "scDRS_run_v8_253_ASC_Pfdr_norm_score",
           "GWAS_sex", "RVAS_sex_fdr001",
           "ASC_Pfdr_253genes_GWAS_vs_fdr001"),
  
  make_cfg(GWAS_Matoba_200, RVAS_fdr001,
           "scDRS_run_GWAS_v3_200_ASD_Matoba_2020_norm_score",
           "scDRS_run_v8_253_ASC_Pfdr_norm_score",
           "GWAS_sex", "RVAS_sex_fdr001",
           "ASC_Pfdr_253genes_200GWAS_vs_fdr001"),
  
  make_cfg(GWAS_Matoba_100, RVAS_fdr01,
           "scDRS_run_GWAS_v3_100_ASD_Matoba_2020_norm_score",
           "scDRS_run_v8_416_ASC_Pfdr_norm_score",
           "GWAS_sex", "RVAS_sex_fdr01",
           "ASC_Pfdr_416genes_GWAS_vs_fdr01"),
  
  make_cfg(GWAS_Matoba_100, RVAS_fdr05,
           "scDRS_run_GWAS_v3_100_ASD_Matoba_2020_norm_score",
           "scDRS_run_v8_696_ASC_Pfdr_norm_score",
           "GWAS_sex", "RVAS_sex_fdr05",
           "ASC_Pfdr_696genes_GWAS_vs_fdr05"),
  
  make_cfg(GWAS_Matoba_100, RVAS_fdr1,
           "scDRS_run_GWAS_v3_100_ASD_Matoba_2020_norm_score",
           "scDRS_run_v8_951_ASC_Pfdr_norm_score",
           "GWAS_sex", "RVAS_sex_fdr1",
           "ASC_Pfdr_951genes_GWAS_vs_fdr1"),

)

# ── Run ────────────────────────────────────────────────────────────────────────
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
for (cfg in comparisons)        run_comparison(cfg, resultdir)