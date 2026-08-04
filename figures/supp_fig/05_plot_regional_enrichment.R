# =============================================================================
# Burden Analysis — O/E ratio with updated gene sets
# Gene sets: GZ=26, CP=5, TH=44 (derived below — see "Assign genes to regional factors")
# Variants: PTV + Mis2 + Mis1
# Cohorts: ASD_DA_Proband (N=24,839) | ASD_D_Proband (N=13,841) | ASD_Proband | Sibling (N=9,567)
# Pairwise: 24 tests total, Bonferroni 0.05/24 = 0.00278
#   Dim1: region pairs within cohort (12 tests)
#   Dim2: cohort pairs within region (12 tests)
# =============================================================================

library(dplyr)
library(tibble)
library(readxl)
library(data.table)
library(openxlsx)

datadir   <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"

# ── Assign genes to regional factors ─────────────────────────────
# Regional NMF factors: Factor 2 = Germinal Zones, Factor 3 = Cortical Plate, Factor 5 = Thalamus
# (Factor 1 = Subplate and Factor 4 = Caudate & Putamen are not spatial regions used here.)
#
# Assignment rule
#   1. Restrict to genes in the RVAS FDR<0.001 universe (full_results_fdr001.txt).
#   2. A gene loading on exactly one of factors {2,3,5} is assigned to that region ("Unique member").
#   3. A gene loading on >1 of factors {2,3,5} is assigned to whichever has the highest
#      relative loading ("Highest loading").
#   4. If the max loading is tied across >=2 of those factors, the gene is a "tie gene" and
#      excluded entirely (e.g. ADNP, TLE3, WAC).
#   5. A gene not loading on any of factors {2,3,5} is not used ("Not in any factor").
#   6. Of the assigned genes, only those flagged Robust expression == TRUE go into the final
#      burden-test gene sets 
region_factor_map <- c(`2` = "Germinal Zones", `3` = "Cortical Plate", `5` = "Thalamus")
region_key_map     <- c("Germinal Zones" = "Germinal_Zones",
                        "Cortical Plate" = "Cortical_Plate",
                        "Thalamus"       = "Thalamus")

# Builds the full gene x regional-factor table — one row per FDR<0.001 gene, with
# Yes/No + loading for each of factors 2/3/5, the resolved assignment, and the
# robust-expression flag. 
build_gene_factor_mapping <- function(enriched_path, fdr001_path) {
  fdr001_genes <- fread(fdr001_path)$gene
  enriched     <- fread(enriched_path)
  spatial      <- enriched %>% filter(factor %in% as.integer(names(region_factor_map)))
  
  assign_one_gene <- function(g) {
    rows <- spatial %>% filter(gene == g)
    if (nrow(rows) == 0) return(list(region = "Not in any factor", reason = NA_character_))
    if (nrow(rows) == 1) return(list(region = region_factor_map[[as.character(rows$factor)]],
                                     reason = "Unique member"))
    max_loading <- max(rows$relative_loading)
    winners     <- rows %>% filter(relative_loading == max_loading)
    if (nrow(winners) > 1) return(list(region = "Excluded (tie)",
                                       reason = "Equal loading ≥2 factors — removed"))
    list(region = region_factor_map[[as.character(winners$factor)]],
         reason = paste0("Highest loading (", winners$relative_loading, ")"))
  }
  assignments <- lapply(fdr001_genes, assign_one_gene)
  
  wide_region <- function(fnum) {
    d <- spatial %>% filter(factor == fnum) %>% select(gene, relative_loading)
    tibble(gene = fdr001_genes) %>%
      left_join(d, by = "gene") %>%
      mutate(yn = if_else(is.na(relative_loading), "No", "Yes"))
  }
  gz <- wide_region(2); cp <- wide_region(3); th <- wide_region(5)
  robust_lookup <- enriched %>% distinct(gene, robust_expression) %>% deframe()
  
  tibble(
    Gene = fdr001_genes,
    `Factor 2\n(GZ)?` = gz$yn,                       `Loading\n(GZ)` = gz$relative_loading,
    `Factor 3\n(CP)?` = cp$yn,                       `Loading\n(CP)` = cp$relative_loading,
    `Factor 5\n(TH)?` = th$yn,                       `Loading\n(TH)` = th$relative_loading,
    `Assigned\nRegion`   = vapply(assignments, `[[`, character(1), "region"),
    `Assignment\nReason` = vapply(assignments, `[[`, character(1), "reason"),
    `Robust\nExpression` = if_else(robust_lookup[Gene], "Yes", "No")
  )
}

gene_factor_mapping <- build_gene_factor_mapping(
  enriched_path = paste0(resultdir, "region_enrichment/enriched_genes_per_factor.txt"),
  fdr001_path   = paste0(datadir, "RVAS_result/full_results_fdr001.txt")
)

gene_sets <- gene_factor_mapping %>%
  filter(`Assigned\nRegion` %in% names(region_key_map), `Robust\nExpression` == "Yes") %>%
  mutate(region_key = unname(region_key_map[`Assigned\nRegion`])) %>%
  group_by(region_key) %>%
  summarise(genes = list(sort(Gene)), .groups = "drop") %>%
  { setNames(.$genes, .$region_key) }
message("Gene set sizes — ", paste(names(gene_sets), lengths(gene_sets), sep="=", collapse=", "))

N <- list(ASD_DA_Proband=24839L, ASD_D_Proband=13841L, ASD_Proband=38680L, Sibling=9567L)


# ── Load data ─────────────────────────────────────────────────────
load_data <- function(nodd_path, ASD_D_Proband_path, full_path, kaplanis_path) {
  nodd <- read_excel(nodd_path) %>%
    mutate(count_ASD_DA_Proband = PTV_Proband + Mis2_Proband + Mis1_Proband,
           count_sib  = PTV_Sibling + Mis2_Sibling + Mis1_Sibling)
  ASD_D_Proband <- read_excel(ASD_D_Proband_path) %>%
    mutate(count_ASD_D_Proband = PTV_Proband + Mis2_Proband + Mis1_Proband)
  # Kaplanis has sex-split columns; sum Male+Female for each variant class
  kaplanis <- fread(kaplanis_path) %>%
    mutate(count_ASD_D_Proband_kap = PTV_Proband_Male + PTV_Proband_Female +
             Mis2_Proband_Male + Mis2_Proband_Female +
             Mis1_Proband_Male + Mis1_Proband_Female)
  full <- fread(full_path) %>%
    rename(Gene = gene) %>%
    mutate(mu_total = mu.lof + mu.mis2 + mu.mis1)
  
  nodd %>%
    select(Gene, count_ASD_DA_Proband, count_sib) %>%
    left_join(ASD_D_Proband %>% select(Gene, count_ASD_D_Proband), by="Gene") %>%
    left_join(kaplanis %>% select(Gene, count_ASD_D_Proband_kap), by="Gene") %>%
    left_join(full %>% select(Gene, mu_total, mu.lof, mu.mis2, mu.mis1), by="Gene") %>%
    mutate(
      count_ASD_D_Proband = count_ASD_D_Proband + coalesce(count_ASD_D_Proband_kap, 0L),
      count_ASD_Proband = count_ASD_DA_Proband + count_ASD_D_Proband
    ) %>%
    select(-count_ASD_D_Proband_kap)
}

# ── O/E with Garwood CI ───────────────────────────────────────────
oe_garwood <- function(obs, exp) {
  oe   <- obs / exp
  ci_l <- if (obs > 0) qchisq(0.025, 2*obs) / (2*exp) else 0
  ci_h <- qchisq(0.975, 2*(obs+1)) / (2*exp)
  list(oe=oe, ci_l=ci_l, ci_h=ci_h)
}

# ── Pairwise conditional binomial ────────────────────────────────
pairwise_test <- function(obs1, exp1, obs2, exp2) {
  obs1 <- as.numeric(obs1)[1]; exp1 <- as.numeric(exp1)[1]
  obs2 <- as.numeric(obs2)[1]; exp2 <- as.numeric(exp2)[1]
  
  if (is.na(obs1) || is.na(obs2) || obs1 == 0 || obs2 == 0)
    return(list(rel=NA, cl=NA, ch=NA, p=NA))
  
  rel <- (obs1/exp1) / (obs2/exp2)
  p0  <- exp1 / (exp1 + exp2)
  p   <- binom.test(obs1, obs1+obs2, p0, alternative="two.sided")$p.value
  se  <- sqrt(1/obs1 + 1/obs2)
  list(rel=rel, cl=exp(log(rel)-qnorm(0.975)*se),
       ch=exp(log(rel)+qnorm(0.975)*se), p=p)
}

FDR_THRESH <- 0.05
sig_label <- function(p_fdr) {
  case_when(
    is.na(p_fdr)         ~ "NA",
    p_fdr < 0.001        ~ "***",
    p_fdr < 0.01         ~ "**",
    p_fdr < FDR_THRESH   ~ "*",
    TRUE                 ~ "ns"
  )
}
# ── Main analysis ─────────────────────────────────────────────────
run_burden <- function(dat) {
  cohort_cols <- list(
    ASD_DA_Proband     = list(col="count_ASD_DA_Proband",     N=N$ASD_DA_Proband),
    ASD_D_Proband     = list(col="count_ASD_D_Proband",     N=N$ASD_D_Proband),
    ASD_Proband = list(col="count_ASD_Proband", N=N$ASD_Proband),
    Sibling  = list(col="count_sib",      N=N$Sibling)
  )
  
  # O/E per region per cohort
  oe_results <- lapply(names(cohort_cols), function(cohort) {
    col <- cohort_cols[[cohort]]$col
    Nc  <- cohort_cols[[cohort]]$N
    lapply(names(gene_sets), function(region) {
      genes <- gene_sets[[region]]
      m     <- dat %>% filter(Gene %in% genes)
      obs   <- sum(m[[col]], na.rm=TRUE)
      exp   <- sum(m$mu_total, na.rm=TRUE) * 2 * Nc
      ci    <- oe_garwood(obs, exp)
      data.frame(cohort=cohort, region=region, n_genes=length(genes),
                 N=Nc, obs=obs, exp=round(exp,4),
                 oe=ci$oe, ci_l=ci$ci_l, ci_h=ci$ci_h)
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  stopifnot(nrow(oe_results) == length(cohort_cols) * length(gene_sets))
  
  # Dim1: region pairs within cohort
  region_pairs <- list(
    c("Cortical_Plate","Germinal_Zones"),
    c("Thalamus","Germinal_Zones"),
    c("Cortical_Plate","Thalamus")
  )
  dim1 <- lapply(names(cohort_cols), function(cohort) {
    lapply(region_pairs, function(pair) {
      ra <- pair[1]; rb <- pair[2]
      sa <- oe_results %>% filter(cohort==!!cohort, region==ra)
      sb <- oe_results %>% filter(cohort==!!cohort, region==rb)
      pw <- pairwise_test(sa$obs, sa$exp, sb$obs, sb$exp)
      data.frame(cohort=cohort,
                 comparison=paste(gsub("_"," ",ra),"vs",gsub("_"," ",rb)),
                 rel_rr=pw$rel, ci_l=pw$cl, ci_h=pw$ch,
                 p=pw$p, sig=sig_label(pw$p))
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  # Dim2: cohort pairs within region
  group_pairs <- list(c("ASD_Proband","Sibling"), c("ASD_D_Proband","ASD_DA_Proband"), c("ASD_D_Proband","Sibling"), c("ASD_DA_Proband","Sibling"))
  dim2 <- lapply(names(gene_sets), function(reg) {  
    lapply(group_pairs, function(pair) {
      ga <- pair[1]; gb <- pair[2]
      sa <- oe_results %>% filter(cohort==ga, region==reg) 
      sb <- oe_results %>% filter(cohort==gb, region==reg)  
      pw <- pairwise_test(sa$obs, sa$exp, sb$obs, sb$exp)
      data.frame(region=gsub("_"," ",reg),                 
                 comparison=paste(ga,"vs",gb),
                 rel_rr=pw$rel, ci_l=pw$cl, ci_h=pw$ch,
                 p=pw$p, sig=sig_label(pw$p))
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  all_p <- c(dim1$p, dim2$p)
  fdr   <- p.adjust(all_p, method = "BH")
  dim1$p_fdr <- fdr[seq_len(nrow(dim1))]
  dim2$p_fdr <- fdr[seq_len(nrow(dim2)) + nrow(dim1)]
  dim1 <- dim1 %>% mutate(sig = sig_label(p_fdr))
  dim2 <- dim2 %>% mutate(sig = sig_label(p_fdr))
  
  list(oe=oe_results, dim1=dim1, dim2=dim2)
}

dat <- load_data(paste0(datadir, "RVAS_result/counts_asd_noddid.xlsx"),
                 paste0(datadir, "RVAS_result/counts_asd_ddid.xlsx"),
                 paste0(datadir, "RVAS_result/full_results_wcounts_2025-10-08.txt"),
                 paste0(datadir, "RVAS_result/Kaplanis_DDD_ASD_counts_by_gene_new_mis_cats_2025-04-09.txt"))
res <- run_burden(dat)
print(res$oe)
print(res$dim1)
print(res$dim2)

fwrite(res$oe, file = paste0(resultdir, "region_enrichment/Obeserved_expected_ratio.csv"),sep = "\t")
fwrite(res$dim1, file = paste0(resultdir, "region_enrichment/Region_comparison.csv"),sep = "\t")
fwrite(res$dim2, file = paste0(resultdir, "region_enrichment/Cohort_comparison.csv"),sep = "\t")

# ── Write manuscript supplementary table ───────────────
# Sheet A = gene_factor_mapping, B = O/E ratios, C = cohort-pair comparisons (per region),
# D = region-pair comparisons (per cohort). 
table_dir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/output/manuscript/tables/"

wb <- createWorkbook()
sheets <- list(
  "A. gene_factor_mapping" = gene_factor_mapping,
  "B. O_E_ratio"           = res$oe,
  "C. Cohort_comparison"   = res$dim2,
  "D. Region_comparison"   = res$dim1
)
header_style <- createStyle(textDecoration = "bold", wrapText = TRUE, valign = "center")
for (sheet_name in names(sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sheets[[sheet_name]])
  addStyle(wb, sheet_name, header_style, rows = 1, cols = seq_len(ncol(sheets[[sheet_name]])))
  setRowHeights(wb, sheet_name, rows = 1, heights = 30)
  freezePane(wb, sheet_name, firstRow = TRUE)
}
saveWorkbook(wb, file = paste0(table_dir, "TableS15_spatial_enrichment.xlsx"), overwrite = TRUE)

# ── plot ──────────────────────────────────────────────────────────
library(ggplot2); library(patchwork); library(dplyr)
COL_ASD_DA_Proband <- "#4472C4"; COL_ASD_D_Proband <- "#C0392B"; COL_ALL <- "#2E7D32"; COL_SIB <-"#7D3C98"
base_theme <- theme_classic(base_size = 6) +
  theme(axis.line=element_line(colour="black",linewidth=0.4),
        axis.ticks=element_line(colour="black",linewidth=0.4),
        axis.text=element_text(colour="black",size=6),
        axis.title=element_text(size=6),
        plot.title=element_text(size=8,hjust=0.5,face = "bold"),
        legend.position="none", panel.grid=element_blank())

region_levels <- c("Thalamus\n(n=44)","Cortical plate\n(n=5)","Germinal zones\n(n=26)")
cohort_levels <- c("Sibling","ASD_Proband","ASD_D_Proband","ASD_DA_Proband")

# ── Panel 1: total observed DNMs ─────────────────────────────────
counts_df <- res$oe %>%
  filter(cohort %in% c("ASD_D_Proband","ASD_DA_Proband","ASD_Proband","Sibling")) %>%
  select(region, cohort, n_genes, obs) %>%
  mutate(
    region = recode(region,
                    "Germinal_Zones" = "Germinal zones\n(n=26)",
                    "Cortical_Plate" = "Cortical plate\n(n=5)",
                    "Thalamus"       = "Thalamus\n(n=44)"),
    region = factor(region, levels = c("Thalamus\n(n=44)",
                                       "Cortical plate\n(n=5)",
                                       "Germinal zones\n(n=26)")),
    cohort = factor(cohort, levels = c("ASD_D_Proband","ASD_DA_Proband","ASD_Proband","Sibling"))
  )

# ── Panel 2: avg DNMs per gene per sample (×10⁻⁴) ─────────────────
avg_df <- res$oe %>%
  filter(cohort %in% c("ASD_D_Proband","ASD_DA_Proband","ASD_Proband","Sibling")) %>%
  mutate(
    avg    = obs / (n_genes * N) * 1e4,
    avg_cl = qchisq(0.025, 2 * obs)     / 2 / (n_genes * N) * 1e4,
    avg_ch = qchisq(0.975, 2 * (obs+1)) / 2 / (n_genes * N) * 1e4,
    region  = recode(region,
                     "Germinal_Zones" = "Germinal zones\n(n=26)",
                     "Cortical_Plate" = "Cortical plate\n(n=5)",
                     "Thalamus"       = "Thalamus\n(n=44)"),
    region  = factor(region, levels = c("Thalamus\n(n=44)",
                                        "Cortical plate\n(n=5)",
                                        "Germinal zones\n(n=26)")),
    cohort  = factor(cohort, levels = c("ASD_D_Proband","ASD_DA_Proband","ASD_Proband","Sibling"))
  ) %>%
  select(region, cohort, n_genes, N, obs, exp, avg, avg_cl, avg_ch)

dodge_width <- 0.80
cohort_levels <- c("ASD_D_Proband","ASD_DA_Proband","ASD_Proband","Sibling")  
n_cohorts <- length(cohort_levels)

avg_df <- avg_df %>%
  mutate(
    cohort_idx = as.integer(factor(cohort, levels = cohort_levels)),
    y_pos = as.numeric(region) +
      (cohort_idx - (n_cohorts + 1) / 2) * (dodge_width / n_cohorts)
  )

# ── Build panels ─────────────────────────────────────────────────
p1 <- ggplot(counts_df,aes(x=obs,y=region,fill=cohort)) +
  geom_col(position=position_dodge(width=0.80),width=0.72) +
  scale_fill_manual(values=c("ASD_DA_Proband"=COL_ASD_DA_Proband,"ASD_D_Proband"=COL_ASD_D_Proband,"ASD_Proband"=COL_ALL, "Sibling"=COL_SIB)) +
  scale_x_continuous(breaks=seq(0,800,200),expand=expansion(mult=c(0,0.05))) +
  labs(title="Total DNM count",x="# DNMs in each gene set",y=NULL) +
  base_theme + theme(axis.text.y=element_text(hjust=0,lineheight=0.9))

p2 <- ggplot(avg_df) +
  geom_col(aes(x=avg,y=region,fill=cohort),
           position=position_dodge(width=0.80), width=0.72) +
  geom_errorbarh(aes(xmin=avg_cl,xmax=avg_ch,y=y_pos),
                 colour="grey25",width=0.05,linewidth=0.3) +
  scale_fill_manual(values=c("ASD_DA_Proband"=COL_ASD_DA_Proband,"ASD_D_Proband"=COL_ASD_D_Proband,"ASD_Proband"=COL_ALL,"Sibling"=COL_SIB)) +
  scale_x_continuous(breaks=seq(0,10,2),expand=expansion(mult=c(0,0.08))) +
  labs(title="Avg. DNMs per gene, per sample",
       x=expression("Avg. DNMs per gene, per sample ("%.%10^{-4}*")"),y=NULL) +
  base_theme + theme(axis.text.y=element_text(hjust=0,lineheight=0.9))


legend_data <- data.frame(
  x=c(1, 3, 4.8, 7), y=1,
  cohort=factor(c("Sibling","ASD_Proband","ASD_DA_Proband","ASD_D_Proband"),
                levels=c("Sibling","ASD_Proband","ASD_DA_Proband","ASD_D_Proband")),
  label=c("Sibling (n=9,567)","Proband (n=38,680)",
          "ASD_DA_Proband (n=24,839)","ASD_D_Proband (n=13,841)"))

legend_plot <- ggplot(legend_data, aes(x, y, fill=cohort)) +
  geom_tile(width=0.28, height=0.25) +
  geom_text(aes(x=x+0.22, label=label),
            hjust=0, vjust=0.5, size=2) +
  scale_fill_manual(values=c("ASD_DA_Proband"=COL_ASD_DA_Proband,"ASD_D_Proband"=COL_ASD_D_Proband,
                             "ASD_Proband"=COL_ALL,"Sibling"=COL_SIB)) +
  xlim(0.6, 10) + ylim(0.7, 1.3) +
  theme_void() + theme(legend.position="none")

fig <- (legend_plot/(p1|p2)) + plot_layout(heights=c(0.09,1))
ggsave(
  filename = paste0(resultdir, "region_enrichment/DMN_count_comparison.pdf"),
  plot   = fig,
  width  = 8,
  height = 4,
  units  = "in"
)


library(ggplot2)
library(dplyr)
df_oe <- res$oe %>%
  mutate(region = recode(region,
                         "Germinal_Zones" = "GZ",
                         "Cortical_Plate" = "CP",
                         "Thalamus" = "THL"),
         region = factor(region,
                         levels = c("GZ", "CP", "THL")),
         cohort = factor(cohort,
                         levels = c("ASD_D_Proband","ASD_DA_Proband","ASD_Proband","Sibling")
         ))
sig_df <- res$dim2 %>%
  mutate(region = recode(region,
                         "Germinal Zones" = "GZ",
                         "Cortical Plate" = "CP",
                         "Thalamus" = "THL"),
         region = factor(region,
                         levels = c("GZ", "CP", "THL"))) %>%
  subset(grepl("vs Sibling", comparison), c("comparison","p_fdr","region"))  %>%
  mutate(
    cohort = factor(case_when(
      grepl("^ASD_D_Proband",   comparison) ~ "ASD_D_Proband",
      grepl("^ASD_DA_Proband", comparison) ~ "ASD_DA_Proband",
      grepl("^Proband",        comparison) ~ "ASD_Proband"
    ), levels = c("ASD_D_Proband", "ASD_DA_Proband", "ASD_Proband", "Sibling")),
    sig_label = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE          ~ "ns"
    )) %>%
  left_join(
    df_oe %>% select(region, cohort, ci_h),
    by = c("region", "cohort")
  )

df_pair <- res$dim1 %>%
  filter(cohort %in% c("ASD_D_Proband", "ASD_DA_Proband")) %>%
  mutate(
    sig_label = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    cohort = factor(
      cohort,
      levels = c("ASD_D_Proband", "ASD_DA_Proband")
    ),
    comparison = recode(
      comparison,
      "Thalamus vs Germinal Zones"       = "THL vs GZ",
      "Cortical Plate vs Germinal Zones" = "CP vs GZ",
      "Cortical Plate vs Thalamus"       = "CP vs THL"
    ),
    comparison = factor(
      comparison,
      levels = c(
        "THL vs GZ",
        "CP vs GZ",
        "CP vs THL"
      )
    )
  )

p_oe <- ggplot(df_oe, aes(x = cohort, y = oe, color = cohort)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_h),
                width = 0.2) +
  geom_text(data = sig_df,
            aes(x = cohort, y = ci_h + 1, label = sig_label),
            size = 3, color = "black") +
  scale_color_manual(values = c(
    "ASD_D_Proband"    = "#C0392B",
    "ASD_DA_Proband"  = "#4472C4",
    "ASD_Proband" = "#2E7D32",
    "Sibling" = "#7D3C98"
  )) +
  facet_wrap(~ region, nrow = 1) +
  labs(x = NULL,
       y = "O/E ratio",
       color = NULL) +
  theme_bw(base_size = 6) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 6, face = "bold"),
    axis.text.x      = element_blank(),
    axis.text.y     = element_text(size = 6),
    legend.text     = element_text(size = 6),
    axis.ticks.x     = element_blank(),
    panel.grid      = element_blank()
  )

p_pair <- ggplot(df_pair,
                 aes(x = rel_rr, y = comparison, color = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbarh(aes(xmin = ci_l, xmax = ci_h),
                 height = 0.2,
                 position = position_dodge(0.5)) +
  geom_text(aes(x = ci_h + 0.05, label = sig_label),
            position = position_dodge(0.5),
            hjust = 0, size = 3, color = "black") +
  scale_color_manual(values = c(
    "ASD_D_Proband"    = "#C0392B",
    "ASD_DA_Proband"  = "#4472C4"
  )) +
  labs(x = "Relative Risk",
       y = NULL,
       color = NULL) +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6,angle = 90,vjust = 0.5,hjust = 0.5),
    axis.text.x     = element_text(size = 6),
    axis.title.x    = element_text(size = 6),
    legend.text     = element_text(size = 6)
  )


plot <- p_oe / p_pair 

# ── Combine DMN-count figure and O/E-ratio figure into a single row ─
final_plot <- (wrap_elements(full = fig) | plot) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
ggsave(
  filename = paste0(resultdir, "region_enrichment/DMN_count_and_OE_ratio.pdf"),
  plot   = final_plot,
  width  = 7.5,
  height = 4,
  units  = "in"
)

# =============================================================================
# Spatial factor (GZ/CP/TH) vs. GO topic (SYN/GR/MORPH): gene counts + Fisher's
# exact test.
# Plot: grouped bar chart of gene counts per topic, stratified by region
#   (GZ/CP/TH all shown), split into robust vs non-robust spatial-factor expression.
# Stats: Fisher's exact test (GZ vs TH only, per topic; CP excluded, underpowered),
#   run separately for robust vs non-robust genes, reported as text only.
# Input: TableS15 (gene -> spatial region), TableS9 (gene -> GO topic)
# =============================================================================

library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)

table_dir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/output/manuscript/tables/"

region_df <- read_excel(paste0(table_dir, "TableS15_spatial_enrichment.xlsx"),
                        sheet = "A. gene_factor_mapping") %>%
  rename(gene_name = Gene, region = `Assigned\nRegion`, robust = `Robust\nExpression`) %>%
  select(gene_name, region, robust)

topic_df <- read_excel(paste0(table_dir, "TableS9_assign-genes-to-GO-topic-subcluster.xlsx"),
                       sheet = "A. GO topic") %>%
  select(gene_name, GO_topic)

dat_all <- topic_df %>%
  left_join(region_df, by = "gene_name") %>%
  filter(region %in% c("Germinal Zones", "Cortical Plate", "Thalamus"))

# GZ vs TH only for the Fisher test (Cortical Plate has n=7, underpowered)
dat <- dat_all %>% filter(region %in% c("Germinal Zones", "Thalamus"))

region_levels <- c("Germinal Zones", "Cortical Plate", "Thalamus")
topic_levels  <- c("GR", "MORPH", "SYN")

# ── Fisher's exact test (Germinal Zones vs Thalamus, per topic) ────────────
# Run separately for robust vs non-robust spatial-factor expression: only the robust
# genes go into the region burden test, so that's the primary result. The non-robust
# version is reported for transparency.
run_fisher_gz_th <- function(data) {
  topics <- unique(data$GO_topic)
  lapply(topics, function(top) {
    a  <- sum(data$region == "Germinal Zones" & data$GO_topic == top)
    b  <- sum(data$region == "Germinal Zones" & data$GO_topic != top)
    cc <- sum(data$region == "Thalamus" & data$GO_topic == top)
    d  <- sum(data$region == "Thalamus" & data$GO_topic != top)
    ft <- fisher.test(matrix(c(a, cc, b, d), nrow = 2))
    data.frame(topic = top,
               n_GZ = a, total_GZ = a + b,
               n_TH = cc, total_TH = cc + d,
               or_GZ_vs_TH = unname(ft$estimate),
               ci_l = ft$conf.int[1], ci_h = ft$conf.int[2],
               p = ft$p.value)
  }) %>% bind_rows() %>%
    mutate(p_fdr = p.adjust(p, method = "BH"),
           sig = case_when(
             p_fdr < 0.001 ~ "***",
             p_fdr < 0.01  ~ "**",
             p_fdr < 0.05  ~ "*",
             TRUE          ~ "ns")) %>%
    arrange(match(topic, topic_levels))
}

fisher_or_robust   <- run_fisher_gz_th(dat %>% filter(robust == "Yes"))

# ── Text report ──────────────────────────────────────────────────────────
report_fisher <- function(fisher_or, label) {
  cat(sprintf("Fisher's exact test: Germinal Zones vs Thalamus, per GO topic (%s)\n", label))
  cat(strrep("-", 60), "\n")
  for (i in seq_len(nrow(fisher_or))) {
    r <- fisher_or[i, ]
    cat(sprintf(
      "%-6s  GZ: %2d/%2d (%.0f%%)   TH: %2d/%2d (%.0f%%)   OR(GZ vs TH) = %s   95%% CI [%.2f, %s]   p = %.3g (FDR p = %.3g, %s)\n",
      r$topic, r$n_GZ, r$total_GZ, 100 * r$n_GZ / r$total_GZ,
      r$n_TH, r$total_TH, 100 * r$n_TH / r$total_TH,
      ifelse(is.infinite(r$or_GZ_vs_TH), "Inf", sprintf("%.2f", r$or_GZ_vs_TH)),
      r$ci_l, ifelse(is.infinite(r$ci_h), "Inf", sprintf("%.2f", r$ci_h)),
      r$p, r$p_fdr, r$sig
    ))
  }
  cat("\n")
}

report_fisher(fisher_or_robust,   "robust spatial-factor expressed only")
report_fisher(fisher_or_unrobust, "non-robust spatial-factor expressed only")

# ── Plot: gene counts per spatial factor, stratified by GO topic ───────────
# Split by robust-expression status: "Yes" is the set that actually goes into the
# region burden test; "No" is assigned-but-filtered-out genes (see TableS11 sheet A).
count_df <- dat_all %>%
  mutate(region   = factor(region, levels = region_levels),
         GO_topic = factor(GO_topic, levels = topic_levels),
         robust   = factor(robust, levels = c("Yes", "No"))) %>%
  count(GO_topic, region, robust, .drop = FALSE)

base_theme <- theme_classic(base_size = 6) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 6),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "top",
        panel.grid = element_blank())

make_count_plot <- function(robust_flag, subtitle) {
  ggplot(count_df %>% filter(robust == robust_flag),
         aes(x = region, y = n, fill = GO_topic)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(aes(label = n), position = position_dodge(width = 0.75),
              vjust = -0.4, size = 2.3) +
    scale_fill_manual(values = c(GR = "#e98024", MORPH = "#157e88", SYN = "#4b2c77")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Number of genes", title = subtitle) +
    base_theme + 
    coord_flip()
}

p_robust   <- make_count_plot("Yes", "Robust spatial-factor expressed\n(used in burden test)")
p_unrobust <- make_count_plot("No",  "Non-robust spatial-factor expressed\n(excluded)")

p <- p_robust + p_unrobust +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

print(p)

ggsave(
  filename = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/region_enrichment/spatial_topic_gene_counts.pdf",
  plot = p, width = 6, height = 3, units = "in")