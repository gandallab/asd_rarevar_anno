library(data.table)
library(dplyr)
library(ggplot2)

# ══════════════════════════════════════════════════════════════════════════════
# Per-DONOR group-level scDRS Z-score, PFC vs V1, pooling all excitatory (or
# inhibitory) neurons across cell types and age groups. Group key = (EN/IN) x donor x
# region, restricted to the 11-donor paired cohort (subclass_donor_region section of
# 01_2_prepare_downstream_metadata.ipynb, run by 05_run_subclass_donor_region_analysis.sh).
# ══════════════════════════════════════════════════════════════════════════════

resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
gene_number <- 253
trait <- "ASC_Pfdr"

group_dir <- file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis",
                       paste0(gene_number, "_genes"), "subclass_donor_region_paired")

dt <- fread(file.path(group_dir, paste0(trait, ".scdrs_group.subclass_donor_region")), sep = "\t")
dt[, c("subclass_group", "donor", "region") := tstrsplit(group, "_", fixed = TRUE)]
dt <- dt[region %in% c("PFC", "V1")]

wide <- dcast(dt, subclass_group + donor ~ region, value.var = c("assoc_mcz", "n_cell"))
setnames(wide, c("assoc_mcz_PFC", "assoc_mcz_V1"), c("z_PFC", "z_V1"))
wide <- wide[!is.na(z_PFC) & !is.na(z_V1)]
wide[, z_diff_V1_minus_PFC := z_V1 - z_PFC]

fwrite(wide, file.path(resultdir, "WangNature/slingshot", "paired_PFC_V1_subclass_donor_Zscore.csv"))
message(nrow(wide), " (subclass_group x donor) rows with both PFC and V1 Z-scores.")

# Paired test across the 11 donors, one test per subclass_group
paired_test <- wide[, {
  tt <- t.test(z_PFC, z_V1, paired = TRUE)
  wt <- wilcox.test(z_PFC, z_V1, paired = TRUE)
  data.table(n_donors = .N,
             mean_diff_V1_minus_PFC = mean(z_diff_V1_minus_PFC),
             t_p = tt$p.value,
             wilcox_p = wt$p.value)
}, by = subclass_group]
fwrite(paired_test, file.path(resultdir, "WangNature/slingshot", "paired_PFC_V1_subclass_donor_Zscore_test.csv"))
print(paired_test)

# plot
plot_dt <- melt(wide, id.vars = c("subclass_group", "donor", "z_diff_V1_minus_PFC"),
                measure.vars = c("z_PFC", "z_V1"), variable.name = "region", value.name = "z")
plot_dt[, region := fifelse(region == "z_PFC", "PFC", "V1")]
plot_dt[, region := factor(region, levels = c("PFC", "V1"))]
plot_dt[, region_x := fifelse(region == "PFC", "PFC\n(anterior)", "V1\n(posterior)")]
plot_dt[, region_x := factor(region_x, levels = c("PFC\n(anterior)", "V1\n(posterior)"))]

mean_dt <- plot_dt[, .(z = mean(z)), by = .(subclass_group, region, region_x)]

label_dt <- copy(paired_test)
label_dt[, t_label := paste0("paired t-test p = ", signif(t_p, 3))]

subclass_colors <- c(EN = "#B2513E", IN = "#3E6BB2")

p <- ggplot(plot_dt, aes(x = region_x, y = z, group = donor, color = subclass_group)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.85) +
  geom_line(data = mean_dt, aes(group = subclass_group), color = "black", linewidth = 1) +
  geom_point(data = mean_dt, aes(group = subclass_group), shape = 18, size = 4.5, color = "black") +
  geom_text(data = label_dt, aes(x = 1, y = Inf, label = t_label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 2, fontface = "italic") +
  scale_color_manual(values = subclass_colors, guide = "none") +
  facet_wrap(~subclass_group, labeller = as_labeller(c(EN = "Excitatory neurons", IN = "Inhibitory neurons"))) +
  labs(x = NULL, y = "scDRS group-level Z-score\n(per donor, pooled across all types/time points)",
       title = "Per-donor scDRS group Z-score, V1 vs PFC (paired-donor cohort)") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    axis.text.x  = element_text(size = 6),
    axis.text.y  = element_text(size = 6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title   = element_text(size = 8, face = "bold"),
    strip.text   = element_text(size = 6, face = "oblique")
  )

ggsave(file.path(resultdir, "WangNature/slingshot", "paired_PFC_V1_subclass_donor_Zscore_dumbbell.pdf"),
       p, width = 5, height = 3, units = "in")
