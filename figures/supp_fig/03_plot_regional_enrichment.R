# =============================================================================
# Burden Analysis — O/E ratio with updated gene sets
# Gene sets: GZ=26, CP=5, TH=44
# Variants: PTV + Mis2 + Mis1
# Cohorts: ASD_DA_Proband (N=24,839) | ASD_D_Proband (N=13,841) | ASD_Proband | Sibling (N=9,567)
# Pairwise: 24 tests total, Bonferroni 0.05/24 = 0.00278
#   Dim1: region pairs within cohort (12 tests)
#   Dim2: cohort pairs within region (12 tests)
# =============================================================================

library(dplyr)
library(readxl)
library(data.table)

# ── Gene sets ─────────────────────────────────────────────────────
gene_sets <- list(
  Germinal_Zones = c("CHD2","CREBBP","EHMT1","GIGYF1","HNRNPD","KMT2A","KMT2C",
                     "MED13","MEIS2","PHF21A","POGZ","RFX3","SETD5","SIN3A","SPEN",
                     "SYNCRIP","TBL1XR1","TCF4","TLK2","TRIP12","UBR5","VEZF1",
                     "WDFY3","XPO1","ZBTB20","ZNF292"),
  Cortical_Plate = c("KDM6B","MEF2C","MYT1L","NR4A2","SATB2"),
  Thalamus = sort(c("AHDC1","ANK2","AP2S1","ARID1B","ASH1L","ASXL3","AUTS2",
                    "BRSK2","CAPRIN1","CSNK1E","CSNK2A1","CTNNB1","CUL3",
                    "DEAF1","DLG4","DNMT3A","DYNC1H1","DYRK1A","FBXO11",
                    "GABBR2","GRIA2","GRIN2B","IRF2BPL","KCNMA1","KCNQ3",
                    "NAA15","NCKAP1","NF1","NRXN1","PACS1","PPP3CA","PRR12",
                    "PSMD12","PTEN","SCN2A","SHANK2","STXBP1","SYNGAP1",
                    "SYT1","TAOK1","TCF7L2","YWHAG","ZMYM2","ZMYND8"))
)

N <- list(ASD_DA_Proband=24839L, ASD_D_Proband=13841L, ASD_Proband=38680L, Sibling=9567L)


# ── Load data ─────────────────────────────────────────────────────
# Adjust paths as needed
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

datadir = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/"
resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
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

# ── plot ──────────────────────────────────────────────────────────
library(ggplot2); library(patchwork); library(dplyr)
COL_ASD_DA_Proband <- "#4472C4"; COL_ASD_D_Proband <- "#C0392B"; COL_ALL <- "#2E7D32"; COL_SIB <-"#7D3C98"
base_theme <- theme_classic(base_size = 7) +
  theme(axis.line=element_line(colour="black",linewidth=0.4),
        axis.ticks=element_line(colour="black",linewidth=0.4),
        axis.text=element_text(colour="black",size=7),
        axis.title=element_text(size=7),
        plot.title=element_text(size=10,hjust=0.5),
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

plot <- p_oe + p_pair +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
ggsave(
  filename = paste0(resultdir, "region_enrichment/Obeserved_expected_ratio.pdf"),
  plot   = plot,
  width  = 3.5,
  height = 2.5,
  units  = "in"
)
