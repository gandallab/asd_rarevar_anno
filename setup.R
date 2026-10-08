#===============================================================================
# Load libraries and functions, set plot theme & colors
# Last updated: 2026-07-31
#===============================================================================

# --- Libraries --- #

library(tidyverse)
library(tidymodels)
library(tidytext)
library(tidygraph)
library(stm)
library(vip)
library(janitor)
library(readxl)
library(writexl)
library(rtracklayer)

library(ggside)
library(ggrepel)
library(ggpubr)
library(ggh4x)
library(ggraph)
library(ggfortify)
library(ggExtra)
library(ggtext)
library(ggbeeswarm)
library(igraph)
library(eulerr)
library(gridExtra)
library(paletteer)
library(RColorBrewer)
library(pals)
library(patchwork)
library(ragg)

library(tidyomics)
library(clusterProfiler)
library(gprofiler2)
library(STRINGdb)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(GO.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationHub)
library(ensembldb)

# mixed-effects models + fan/wedge (annular sector) plots
library(lme4)
library(lmerTest)
library(emmeans)
library(ggforce)

library(metafor) # pooled effects estimates


count <- dplyr::count
select <- dplyr::select
filter <- dplyr::filter


# --- Set directories ---#

project_dir <- paste0(here::here(), "/")
data_dir <- paste0(project_dir, "data/")
functions_dir <- paste0(project_dir, "functions/")

out_dir <- paste0(project_dir, "outputs/")
objects_dir <- paste0(out_dir, "objects/")
analysis_dir <- paste0(out_dir, "analysis/")
figures_dir <- paste0(out_dir, "figures/")
tables_dir <- paste0(out_dir, "tables/")


# --- Load functions --- #

# source only .R files (skip subdirectories like archive/ and non-R files)
function_scripts <- list.files(functions_dir, pattern = "\\.R$", full.names = TRUE)
walk(function_scripts, source)


# --- Set plot theme --- #

theme_set(
    theme_bw() +
        theme(
            #text = element_text(family = "Arial"),
            plot.title = element_text(size = 9, hjust = 0.5),
            plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray30"),
            plot.tag = element_text(face = "bold", size = 10),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            strip.text = element_text(size = 8),
            legend.title = element_text(size = 7, hjust = 0.5),
            legend.text = element_text(size = 6)
        )
)


# --- Plot settings (reused across scripts) --- #

# Assign colors to GO topics
topic_colors <- c(
    "SYN" = "#4a2377",
    "GR" = "#ea801c",
    "MORPH" = "#0d7d87"
)

## Fig 2F BrainSpan development:
# stage_levels defined by Kang et al., 2011
# 13 developmental periods following Kang et al. 2011 (Table 1)
# Periods 2+3 share the label "Early fetal"; periods 4+5 share "Early mid-fetal".
# Postnatal ages converted to pcw assuming full-term birth at 40 pcw,
stage_levels <- c(
    "Embryonic",
    "Early fetal",
    "Early mid-fetal",
    "Late mid-fetal",
    "Late fetal",
    "Neonatal and early infancy",
    "Late infancy",
    "Early childhood",
    "Middle and late childhood",
    "Adolescence",
    "Young adulthood",
    "Middle adulthood",
    "Late adulthood"
)


