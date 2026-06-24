# preprocess for RNA seq data
library(tidyverse)
library(readxl)
source("utils.R")

# output dir
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# -- load and clean ------------------------------------------------------------
heat_raw = read_excel("genome_transcriptome_analysis/RNAseq_experiment/results_pilot3_w_KO_PFAM.xlsx") %>%
  mutate(
    significant = !is.na(padj) & padj < 0.05,
    Level_B     = na_if(Level_B, "NA"),
    Level_C     = na_if(Level_C, "NA"),
    Pfam_ID     = na_if(Pfam_ID, "NA"),
    Pfam_desc   = na_if(Pfam_desc, "NA")
  )

# one row per gene with direction label
all_genes = heat_raw %>%
  distinct(Geneid, .keep_all = TRUE) %>%
  mutate(
    has_annot = !is.na(Level_B),
    direction = case_when(
      significant & log2FoldChange > 0 ~ "up",
      significant & log2FoldChange < 0 ~ "down",
      TRUE ~ "ns"
    )
  )

# disease/unclassified level b filt
heat_filt = heat_raw %>% filter(!Level_B %in% EXCLUDE_B)

# sig genes distinct per gene x level_c
sig_by_pathway = heat_filt %>%
  filter(!is.na(padj), padj < 0.05) %>%
  distinct(Geneid, Level_C, .keep_all = TRUE)

# genome-wide null proportion
null_genes    = sig_by_pathway %>% distinct(Geneid, .keep_all = TRUE)
p_up_global   = mean(null_genes$log2FoldChange > 0)
p_down_global = 1 - p_up_global
sig_genes = all_genes %>% filter(!is.na(padj), padj < 0.05)

# outputs
saveRDS(all_genes, file.path(DAT, "all_genes.rds"))
saveRDS(sig_genes,      file.path(DAT, "sig_genes.rds"))
saveRDS(heat_filt,      file.path(DAT, "heat_filt.rds"))
saveRDS(sig_by_pathway, file.path(DAT, "sig_by_pathway.rds"))
saveRDS(list(p_up = p_up_global, p_down = p_down_global),
         file.path(DAT, "null_props.rds"))
