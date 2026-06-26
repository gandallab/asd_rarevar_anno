load_syngo_resources <- function(syngopath) {
  syngo.genes <- as.data.frame(readxl::read_excel(file.path(syngopath, "syngo_genes.xlsx")))
  rownames(syngo.genes) <- syngo.genes$ensembl_id

  syngo.annot <- as.data.frame(readxl::read_excel(file.path(syngopath, "syngo_annotations.xlsx")))
  syngo.annot$ensembl_id <- syngo.genes$ensembl_id[
    match(syngo.annot$hgnc_id, syngo.genes$hgnc_id)
  ]

  syngo.ontol <- as.data.frame(readxl::read_excel(file.path(syngopath, "syngo_ontologies.xlsx")))
  rownames(syngo.ontol) <- syngo.ontol$id
  syngo.ontol[syngo.ontol$domain == "BP", "orig.order"] <- seq_len(sum(syngo.ontol$domain == "BP"))
  syngo.ontol[syngo.ontol$domain == "CC", "orig.order"] <- seq_len(sum(syngo.ontol$domain == "CC"))

  list(
    syngo.genes = syngo.genes,
    syngo.annot = syngo.annot,
    syngo.ontol = syngo.ontol
  )
}

get_hotnet_definitions <- function(syngopath) {
  df_cmu_cc <- readxl::read_excel(file.path(syngopath, "Table S1.xlsx"), sheet = "CC")
  df_cmu_bp <- readxl::read_excel(file.path(syngopath, "Table S1.xlsx"), sheet = "BP")

  cc_hotnet <- list(
    "Active zone" = df_cmu_cc %>% dplyr::filter(orig.order %in% c(10, 11, 12, 15)) %>% dplyr::pull(id),
    "Presynaptic membrane" = df_cmu_cc %>% dplyr::filter(orig.order %in% c(33, 36)) %>% dplyr::pull(id),
    "Postsynaptic specialization" = df_cmu_cc %>% dplyr::filter(orig.order %in% c(57, 59, 61, 62, 63, 64, 67)) %>% dplyr::pull(id)
  )

  bp_hotnet <- list(
    "Presynaptic ion channels" = df_cmu_bp %>% dplyr::filter(orig.order %in% c(7, 8, 9)) %>% dplyr::pull(id),
    "Postsynaptic organization" = df_cmu_bp %>% dplyr::filter(orig.order %in% c(108, 109)) %>% dplyr::pull(id)
  )

  list(cc_hotnet = cc_hotnet, bp_hotnet = bp_hotnet)
}

get_hotnet_genes <- function(hotnet_list, syngo.ontol) {
  purrr::map(hotnet_list, function(ids) {
    unique(unlist(strsplit(syngo.ontol[ids, "hgnc_symbol"], ", ", fixed = TRUE)))
  })
}

get_hotnet_gs_order <- function() {
  c(
    "Active zone",
    "Presynaptic membrane",
    "Presynaptic ion channels",
    "Postsynaptic specialization",
    "Postsynaptic organization"
  )
}

get_aucell_cell_order <- function() {
  c(
    "RG-vRG", "RG-tRG", "RG-oRG",
    "IPC-EN", "EN-Newborn", "EN-IT-Immature",
    "EN-L2_3-IT", "EN-L4-IT", "EN-L5-IT", "EN-L6-IT",
    "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP", "EN-L6-CT", "EN-L6b",
    "IPC-Glia",
    "IN-dLGE-Immature", "IN-CGE-Immature", "IN-MGE-Immature",
    "IN-CGE-VIP", "IN-CGE-SNCG", "IN-Mix-LAMP5",
    "IN-MGE-SST", "IN-MGE-PV",
    "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous",
    "OPC", "Oligodendrocyte-Immature", "Oligodendrocyte",
    "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown"
  )
}

get_hotnet_gs_levels <- function() {
  list(
    "Active zone"                 = c("Low", "Medium", "High"),
    "Presynaptic ion channels"    = c("Low", "Medium", "High"),
    "Presynaptic membrane"        = c("Low", "High"),
    "Postsynaptic specialization" = c("Low", "High"),
    "Postsynaptic organization"   = c("Low", "High")
  )
}

get_aucell_celltype_colors <- function() {
  c(
    "RG-vRG"             = "#525252", "RG-tRG"              = "#969696",
    "RG-oRG"             = "#d9d9d9", "IPC-EN"              = "#ffd700",
    "EN-Newborn"         = "#ffa500", "EN-IT-Immature"      = "#ff7f00",
    "EN-L2_3-IT"         = "#ff4500", "EN-L4-IT"            = "#e31a1c",
    "EN-L5-IT"           = "#bd0026", "EN-L6-IT"            = "#800026",
    "EN-Non-IT-Immature" = "#f768a1", "EN-L5-ET"            = "#c51b8a",
    "EN-L5_6-NP"         = "#7a0177", "EN-L6-CT"            = "#49006a",
    "EN-L6b"             = "#2d0040", "IPC-Glia"            = "#1d91c0",
    "IN-dLGE-Immature"   = "#225ea8", "IN-CGE-Immature"     = "#253494",
    "IN-CGE-VIP"         = "#0d3b8e", "IN-CGE-SNCG"         = "#08306b",
    "IN-Mix-LAMP5"       = "#006d2c", "IN-MGE-Immature"     = "#31a354",
    "IN-MGE-SST"         = "#74c476", "IN-MGE-PV"           = "#bae4b3"
  )
}

get_aucell_lineage_colors <- function() {
  c(
    "RG"    = "#969696", "EN" = "#e31a1c",
    "IN"    = "#1f78b4", "Other" = "#b2df8a"
  )
}
