library(data.table)
library(dplyr)
library(ggplot2)
count_DMN <- fread(file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_group.txt")
count_perproband <- fread(file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_perproband_group.txt")

count_combined <- count_DMN %>% select(Gene, count_new_DMN_ddid, count_new_DMN_noddid, p_hat, Group_phat = Group) %>%
  left_join(count_perproband %>% select(Gene, rate_ddid, rate_noddid,Group_perproband),by = "Gene") %>%
  mutate(
    d_rate = rate_ddid - rate_noddid,
    Group_combined = case_when(
      Group_perproband == "DDID"   & Group_phat == "DDID"   ~ "ASD-NDD",
      Group_perproband == "NoDDID" & Group_phat == "NoDDID" ~ "ASD-p",
      Group_perproband == "DDID"   & Group_phat == "NoDDID" ~ "ASD-mixed",
      Group_perproband == "NoDDID" & Group_phat == "DDID"   ~ "NoDDID(pp) / DDID(phat)"
    )
  )
fwrite(count_combined, file = "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/RVAS_result/full_results_fdr001_PtvMis2Del_combined_annotation.txt", sep = "\t")
ggplot(count_combined, 
       aes(x = d_rate, y = p_hat, color = Group_combined)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = sum(count_DMN$count_new_DMN_ddid) / 
               (sum(count_DMN$count_new_DMN_ddid) + sum(count_DMN$count_new_DMN_noddid)),
             linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c(
      "ASD-NDD"         = "#66C2A5",
      "ASD-p"           = "#FDB462",
      "ASD-mixed"  = "grey"
    )
  ) +
  labs(
    x     = "Difference of mutation rate (DDID - NoDDID)",
    y     = "p_hat (DDID mutation fraction)",
    title = "Mutation rate vs p_hat by gene (PTV+Mis2+DEL)",
    color = "Group"
  ) +
  theme_classic(base_size = 8)

ggsave(
  filename = file.path(
    datadir, "RVAS_result/geneset_comparison.pdf"
  ),
  width  = 5,
  height = 5.2,
  units  = "in"
)

group_counts <- count_combined %>%
  count(Group_combined) %>%
  mutate(label = paste0(Group_combined, "\n(n=", n, ")"))

label_map <- setNames(group_counts$label, group_counts$Group_combined)

ggplot(count_combined %>% mutate(Group_label = label_map[Group_combined]), 
       aes(x = p_hat, fill = Group_label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = sum(count_DMN$count_new_DMN_ddid) / 
               (sum(count_DMN$count_new_DMN_ddid) + sum(count_DMN$count_new_DMN_noddid)),
             linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = setNames(
      c("#B2182B", "#2166AC", "grey", "#92C5DE"),
      label_map[c("ASD-NDD", "ASD-p", "ASD-mixed", "ASD-p-lowburden")]
    )
  ) +
  scale_y_reverse() +
  coord_flip() +
  labs(
    x     = "p_hat (DDID mutation fraction)",
    y     = "n_genes",
    title = "p_hat distribution by group (PTV+Mis2+DEL)",
    fill  = "Group"
  ) +
  theme_classic(base_size = 9)
