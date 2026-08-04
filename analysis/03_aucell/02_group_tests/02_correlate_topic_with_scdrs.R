library(data.table)
library(dplyr)
library(purrr)
library(readxl)

source(file.path("analysis", "_shared", "cluster_labels.R"))
source(file.path("analysis", "_shared", "scdrs_aucell_helpers.R"))

resultdir <- "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/"
aucdir <- paste0(resultdir, "WangNature/AUCell/")
gene_number <- 253
trait <- "ASC_Pfdr"
scdrsdir <- paste0(resultdir, "WangNature/scDRS/run_v8/score_file/", gene_number, "_genes/")
group_dir <- file.path(resultdir, "WangNature/scDRS/run_v8/downstream_analysis", paste0(gene_number, "_genes"))

# ============================================================
# Global definitions
# ============================================================
gene_module <- read_excel("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SynGO/2026.06.30_FINAL-gene-topic-assignment.xlsx")
gene_module %>%
  count(GO_topic, name = "n_genes") %>%
  arrange(desc(n_genes))

all_genesets <- gene_module %>%
  group_by(GO_topic) %>%
  summarise(genes = list(gene_name), .groups = "drop") %>%
  tibble::deframe()

gs_order <- names(all_genesets)
timepoint_order <- get_scdrs_timepoint_order()

# ============================================================
# Load data
# ============================================================
auc_cell <- fread("/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/result/WangNature/AUCell/topic/auc_cell.csv")

meta <- fread(
  "/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/metadata.csv"
) %>%
  rename(cell_id = ID, cell_type = type, donor_id = Ident)

ctrl_cols <- get_scdrs_ctrl_cols()

scdrs_scores_raw <- fread(paste0(scdrsdir, "ASC_Pfdr.full_score.gz")) %>%
  as.data.frame() %>%
  rename(cell_id = V1) %>%
  select(cell_id, norm_score, all_of(ctrl_cols))

en_labels <- get_en_labels()
in_labels <- get_in_labels()
lineage_df <- get_lineage_df_from_mclust(resultdir, en_labels, in_labels)

# ============================================================
# Build cor_input: AUC + scDRS + metadata
# ============================================================
whole_set <- fread(file.path(group_dir, paste0(trait, ".scdrs_group.mclust.csv")), sep = "\t")
sig_groups <- unique(whole_set[assoc_mcp_fdr < 0.05, group])

cor_input <- auc_cell %>%
  left_join(lineage_df, by = "cell_id") %>%
  left_join(meta %>% select(cell_id, Group), by = "cell_id") %>%
  left_join(scdrs_scores_raw, by = "cell_id") %>%
  filter(!is.na(norm_score), !is.na(lineage_label), mclust_group %in% sig_groups) %>%
  mutate(Group = factor(Group, levels = timepoint_order))

ct_order <- intersect(get_scdrs_canonical_order(), unique(cor_input$lineage_label))

write.csv(
  cor_input,
  paste0(aucdir, "scDRS_topic_AUCell_cor_input.csv"),
  row.names = FALSE
)

# ============================================================
# Step 0. Overall MC-corrected Spearman correlation
#         lineage_label x gene_set
# ============================================================
cor_overall_ct <- purrr::map_dfr(gs_order, function(gs_name) {
  message("Overall (cell type): ", gs_name)
  purrr::map_dfr(ct_order, function(ct) {
    dat <- as.data.frame(cor_input) %>% filter(lineage_label == ct)
    if (nrow(dat) < 150) return(NULL)
    if (!gs_name %in% colnames(dat)) return(NULL)

    x <- dat[, gs_name]
    if (sd(x, na.rm = TRUE) == 0) return(NULL)

    r_disease <- cor(dat[["norm_score"]], x, method = "spearman", use = "complete.obs")
    if (is.na(r_disease)) return(NULL)

    r_ctrl <- sapply(ctrl_cols, function(cc) {
      cor(dat[[cc]], x, method = "spearman", use = "complete.obs")
    })

    data.frame(
      lineage_label = ct,
      n_cells = nrow(dat),
      r_disease = r_disease,
      mc_pval = (sum(abs(r_ctrl) >= abs(r_disease)) + 1) / (length(r_ctrl) + 1),
      gene_set = gs_name
    )
  })
}) %>%
  mutate(
    FDR = p.adjust(mc_pval, method = "BH"),
    r_plot = ifelse(mc_pval < 0.05, r_disease, NA),
    lineage_label = factor(lineage_label, levels = ct_order),
    gene_set = factor(gene_set, levels = gs_order)
  )

write.csv(
  cor_overall_ct,
  paste0(aucdir, "scDRS_topic_AUCell_cor_overall_celltype.csv"),
  row.names = FALSE
)

# ============================================================
# Step 1. MC-corrected Spearman correlation
#         lineage_label x Group x gene_set
# ============================================================
cor_results <- purrr::map_dfr(gs_order, function(gs_name) {
  message("Processing (cell type level): ", gs_name)

  purrr::map_dfr(
    split(
      cor_input,
      list(as.character(cor_input$lineage_label), as.character(cor_input$Group))
    ),
    function(dat) {
      if (nrow(dat) < 150) return(NULL)

      r_disease <- cor(
        dat[["norm_score"]],
        dat[[gs_name]],
        method = "spearman",
        use = "complete.obs"
      )
      r_ctrl <- sapply(ctrl_cols, function(cc) {
        cor(dat[[cc]], dat[[gs_name]], method = "spearman", use = "complete.obs")
      })
      mc_pval <- (sum(abs(r_ctrl) >= abs(r_disease)) + 1) / (length(r_ctrl) + 1)

      data.frame(
        lineage_label = unique(dat$lineage_label),
        Group = unique(dat$Group),
        n_cells = nrow(dat),
        r_disease = r_disease,
        mc_pval = mc_pval,
        gene_set = gs_name
      )
    }
  )
}) %>%
  filter(!is.na(r_disease)) %>%
  mutate(
    FDR = p.adjust(mc_pval, method = "BH"),
    r_plot = ifelse(mc_pval < 0.05, r_disease, NA),
    lineage_label = factor(lineage_label, levels = ct_order),
    gene_set = factor(gene_set, levels = gs_order),
    Group = factor(Group, levels = timepoint_order)
  )

write.csv(
  cor_results,
  paste0(aucdir, "scDRS_topic_AUCell_cor_celltype.csv"),
  row.names = FALSE
)
