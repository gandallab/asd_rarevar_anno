get_en_labels <- function() {
  c(
    "1" = "EN-Newborn",
    "2" = "EN-L4-IT",
    "3" = "EN-Newborn_IPC-EN",
    "4" = "EN-Non-IT-Immature",
    "5" = "EN-IT-Immature",
    "6" = "EN-L2_3-IT",
    "7" = "EN-L4-IT-V1",
    "8" = "vRG_oRG",
    "9" = "vRG",
    "10" = "EN-L2_3_4-IT",
    "11" = "EN-L6-IT",
    "12" = "EN-Newborn",
    "13" = "EN-L5_6-NP_EN-L5-ET_EN-Non-IT-Immature",
    "14" = "EN-L5-IT",
    "15" = "EN-L6-CT_EN-Non-IT-Immature_EN-L6b",
    "16" = "EN-Newborn_EN-IT-Immature",
    "17" = "EN-IT-Immature_EN-L2_3-IT",
    "18" = "EN-L4_5-IT",
    "19" = "EN-L6b",
    "21" = "EN-Non-IT-Immature",
    "22" = "IPC-EN",
    "23" = "EN-L6-CT_EN-Non-IT-Immature",
    "24" = "EN-L6b"
  )
}

get_in_labels <- function() {
  c(
    "1" = "IN-MGE-PV",
    "2" = "IN-CGE-VIP",
    "3" = "IN-MGE-Immature_IN-MGE-SST-upper",
    "4" = "vRG",
    "5" = "IN-MGE-SST-upper",
    "6" = "IN-MGE-Immature",
    "7" = "vRG",
    "8" = "vRG",
    "9" = "IN-MGE-SST-deep",
    "10" = "IN-CGE-LAMP5",
    "11" = "IN-MGE-SST-deep",
    "12" = "IN-CGE-Immature_IN-dLGE-Immature",
    "13" = "vRG",
    "14" = "IN-dLGE-Immature",
    "15" = "IN-MGE-PV",
    "16" = "IN-MGE-PV",
    "17" = "oRG_tRG",
    "18" = "IN-CGE-SNCG",
    "19" = "vRG",
    "20" = "Tri-IPC",
    "22" = "IN-CGE-VIP_IN-CGE-Immature"
  )
}

get_lineage_cluster_labels <- function() {
  list(
    EN = get_en_labels(),
    IN = get_in_labels()
  )
}


get_marker_reference <- function() {
  list(
    "RG-vRG"                 = c("HES1", "TFAP2C"),
    "RG-tRG"                 = c("HES1", "TFAP2C", "CRYAB", "FBXO32"),
    "RG-oRG"                 = c("HES1", "TFAP2C", "TNC", "HOPX"),
    "IPC-EN"                 = c("SLC17A6", "EOMES"),
    "EN-Newborn"              = c("SLC17A6", "NRP1", "CUX2"),
    "EN-IT-Immature"         = c("SLC17A6", "GLIS3", "CUX2"),
    "EN-L2_3-IT"             = c("SLC17A6", "CUX2", "NWD2"),
    "EN-L4-IT"               = c("SLC17A6", "CUX2", "IL1RAPL2", "NWD2", "TSHZ2", "RORB"),
    "EN-L5-IT"               = c("SLC17A6", "RORB", "IL1RAPL2", "FEZF2", "NWD2"),
    "EN-L6-IT"               = c("SLC17A6", "IL1RAPL2", "ZNF804B", "NWD2"),
    "EN-Non-IT-Immature"     = c("SLC17A6", "FEZF2"),
    "EN-L5-ET"               = c("SLC17A6", "RORB", "FEZF2", "NWD2"),
    "EN-L5_6-NP"             = c("SLC17A6", "FEZF2", "TSHZ2"),
    "EN-L6-CT"               = c("SLC17A6", "FEZF2", "SYT6"),
    "EN-L6b"                 = c("SLC17A6", "IL1RAPL2", "FEZF2", "FAM160A1"),
    "IN-dLGE-Immature"       = c("GAD1", "MEIS2"),
    "IN-CGE-Immature"        = c("GAD1", "ADARB2"),
    "IN-CGE-VIP"             = c("GAD1", "ADARB2", "VIP", "CCK"),
    "IN-CGE-SNCG"            = c("GAD1", "ADARB2", "CCK"),
    "IN-Mix-LAMP5"           = c("GAD1", "ADARB2", "LAMP5"),
    "IN-MGE-Immature"        = c("GAD1", "NXPH1"),
    "IN-MGE-SST"             = c("GAD1", "NXPH1", "SST"),
    "IN-MGE-PV"              = c("GAD1", "NXPH1", "PVALB"),
    "IPC-Glia"               = c("EGFR", "MAP3K1", "OLIG1"),
    "Astrocyte-Immature"     = c("AQP4", "MOXD1"),
    "Astrocyte-Protoplasmic" = c("AQP4", "GRM3"),
    "Astrocyte-Fibrous"      = c("AQP4", "GFAP"),
    "OPC"                    = c("MAP3K1", "OLIG1", "SOX10"),
    "Oligodendrocyte-Immature" = c("OLIG1", "SOX10", "BCAS1", "MBP"),
    "Oligodendrocyte"        = c("OLIG1", "SOX10", "GRM3", "MBP"),
    "Cajal-Retzius"          = c("RELN"),
    "Microglia"              = c("IRF8"),
    "Vascular"               = c("IGFBP7")
  )
}

get_marker_gene_types <- function(genes, marker_ref = get_marker_reference()) {
  sapply(genes, function(g) {
    types <- names(marker_ref)[sapply(marker_ref, function(gs) g %in% gs)]
    if (length(types) == 0) "?" else paste(types, collapse = " / ")
  })
}