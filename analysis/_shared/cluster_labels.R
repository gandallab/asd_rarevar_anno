get_en_labels <- function() {
  c(
    "1" = "EN-Newborn",
    "2" = "EN-L4-IT",
    "3" = "EN-Newborn",
    "4" = "EN-Non-IT-Immature",
    "5" = "EN-IT-Immature",
    "6" = "EN-L2_3-IT",
    "7" = "EN-L4-IT-V1",
    "8" = "oRG and tRG",
    "9" = "vRG",
    "10" = "EN-L2_3_4-IT",
    "11" = "EN-L6-IT",
    "12" = "EN-Newborn",
    "13" = "EN-L5_6-NP & EN-L5-ET",
    "14" = "EN-L5-IT",
    "15" = "EN-L6-CT",
    "16" = "EN-Newborn",
    "17" = "EN-IT-Immature",
    "18" = "EN-L4_5-IT",
    "19" = "EN-L6b",
    "21" = "EN-Non-IT-Immature",
    "22" = "IPC-EN",
    "23" = "EN-L6-CT",
    "24" = "EN-L6b"
  )
}

get_in_labels <- function() {
  c(
    "1" = "IN-MGE-PV",
    "2" = "IN-CGE-VIP",
    "3" = "IN-MGE-Immature",
    "4" = "vRG",
    "5" = "IN-MGE-SST-upper",
    "6" = "IN-MGE-Immature",
    "7" = "vRG",
    "8" = "vRG",
    "9" = "IN-MGE-SST-deep",
    "10" = "IN-CGE-LAMP5",
    "11" = "IN-MGE-SST-deep",
    "12" = "IN-CGE-Immature",
    "13" = "vRG",
    "14" = "IN-dLGE-Immature",
    "15" = "IN-MGE-PV",
    "16" = "IN-MGE-PV",
    "17" = "oRG and tRG",
    "18" = "IN-CGE-SNCG",
    "19" = "vRG",
    "20" = "Tri-IPC",
    "22" = "IN-CGE-VIP"
  )
}

get_lineage_cluster_labels <- function() {
  list(
    EN = get_en_labels(),
    IN = get_in_labels()
  )
}

