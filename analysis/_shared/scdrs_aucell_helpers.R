get_scdrs_celltype_groups <- function() {
  list(
    "EN-dev" = c("IPC-EN", "EN-Newborn"),
    "EN-IT" = c("EN-L2_3-IT", "EN-L2_3_4-IT", "EN-L4-IT-V1", "EN-L4-IT", "EN-L5-IT"),
    "IN-dev" = c("IN-dLGE-Immature", "IN-CGE-Immature", "IN-MGE-Immature"),
    "IN-CGE" = c("IN-CGE-LAMP5"),
    "IN-MGE" = c("IN-MGE-SST-upper", "IN-MGE-PV")
  )
}

get_scdrs_ct_order <- function(celltype_groups = get_scdrs_celltype_groups()) {
  unname(unlist(celltype_groups))
}

get_scdrs_group_order <- function() {
  c("EN-dev", "EN-IT", "IN-dev", "IN-CGE", "IN-MGE")
}

get_scdrs_group_map <- function(celltype_groups = get_scdrs_celltype_groups()) {
  purrr::imap_dfr(celltype_groups, function(types, grp) {
    data.frame(lineage_label = types, cell_group = grp)
  })
}

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
    file_in %>% dplyr::select(cell_id, lineage_label = label),
    file_en %>% dplyr::select(cell_id, lineage_label = label)
  ) %>%
    dplyr::distinct(cell_id, .keep_all = TRUE)
}
