get_scdrs_timepoint_order <- function() {
  c("First_trimester", "Second_trimester", "Third_trimester", "Infancy", "Adolescence")
}

get_scdrs_timepoint_colors <- function() {
  c(
    "First_trimester" = "#f0f921",
    "Second_trimester" = "#fca636",
    "Third_trimester" = "#e16462",
    "Infancy" = "#b12a90",
    "Adolescence" = "#6a00a8"
  )
}

get_scdrs_timepoint_labels <- function() {
  c(
    "First_trimester" = "T1",
    "Second_trimester" = "T2",
    "Third_trimester" = "T3",
    "Infancy" = "Inf",
    "Adolescence" = "Ado"
  )
}

get_scdrs_ctrl_cols <- function() {
  paste0("ctrl_norm_score_", 0:999)
}

get_scdrs_canonical_order <- function() {
  c(
    "vRG", "vRG_oRG", "oRG_tRG", "vRG_oRG_tRG",
    "IPC-EN", "EN-Newborn_IPC-EN", "EN-Newborn", "EN-Newborn_EN-IT-Immature",
    "EN-IT-Immature", "EN-IT-Immature_EN-L2_3-IT", "EN-L2_3-IT", "EN-L2_3_4-IT",
    "EN-L4-IT", "EN-L4-IT-V1", "EN-L4_5-IT", "EN-L5-IT", "EN-L6-IT",
    "EN-Non-IT-Immature", "EN-L5_6-NP_EN-L5-ET_EN-Non-IT-Immature",
    "EN-L6-CT_EN-Non-IT-Immature", "EN-L6-CT_EN-Non-IT-Immature_EN-L6b", "EN-L6b",
    "Tri-IPC",
    "IN-dLGE-Immature", "IN-CGE-Immature_IN-dLGE-Immature", "IN-CGE-Immature", "IN-MGE-Immature", 
    "IN-MGE-Immature_IN-MGE-SST-upper", "IN-MGE-SST-upper",
    "IN-MGE-SST-deep", "IN-MGE-PV","IN-CGE-LAMP5", "IN-CGE-VIP_IN-CGE-Immature", "IN-CGE-VIP","IN-CGE-SNCG" 
  )
}

get_lineage_df_from_mclust <- function(resultdir, en_labels, in_labels) {
  file_in <- data.table::fread(file.path(resultdir, "WangNature", "slingshot", "IN", "IN_lineage_mclust.csv")) %>%
    as.data.frame() %>%
    dplyr::rename(cell_id = ID) %>%
    dplyr::mutate(label = in_labels[as.character(mclust22)])
  
  file_en <- data.table::fread(file.path(resultdir, "WangNature", "slingshot", "EN", "EN_lineage_mclust.csv")) %>%
    as.data.frame() %>%
    dplyr::rename(cell_id = ID) %>%
    dplyr::mutate(label = en_labels[as.character(mclust25)])
  
  dplyr::bind_rows(
    file_in %>% dplyr::select(cell_id, lineage_label = label, mclust_group = mclust),
    file_en %>% dplyr::select(cell_id, lineage_label = label, mclust_group = mclust)
  ) %>%
    dplyr::distinct(cell_id, .keep_all = TRUE)
}