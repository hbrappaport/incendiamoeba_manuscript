# utils for themes, data, filtering

library(tidyverse)
library(scales)

OUT      = "genome_transcriptome_analysis/RNAseq_experiment/outputs"
DAT      = "genome_transcriptome_analysis/RNAseq_experiment/rds"
COL_UP   = "#C1532D"
COL_DOWN = "#3A6AA1"

SAMPLE_COLS = c("48_A","48_B","48_C","48_D","61_A","61_B","61_C","61_D")

SAMPLE_INFO = tibble(
  sample    = SAMPLE_COLS,
  condition = c("48C","48C","48C","48C","61C","61C","61C","61C")
)

theme_rna = function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size - 3, colour = "grey30"),
      axis.text     = element_text(colour = "black"),
      axis.title    = element_text(colour = "black"),
      legend.text   = element_text(colour = "black"),
      plot.margin   = margin(8, 12, 8, 8)
    )
}

# categories to exclude/filter
EXCLUDE_B = c(
  "09176 Drug resistance: antineoplastic",
  "09175 Drug resistance: antimicrobial",
  "09174 Infectious disease: parasitic",
  "09172 Infectious disease: viral",
  "09171 Infectious disease: bacterial",
  "09167 Endocrine and metabolic disease",
  "09166 Cardiovascular disease",
  "09165 Substance dependence",
  "09164 Neurodegenerative disease",
  "09163 Immune disease",
  "09162 Cancer: specific types",
  "09161 Cancer: overview",
  "09151 Immune system",
  "09152 Endocrine system",
  "09149 Aging",
  "09158 Development and regeneration",
  "09153 Circulatory system",
  "09154 Digestive system",
  "09156 Nervous system",
  "09155 Excretory system",
  "09157 Sensory system",
  "09159 Environmental adaptation",
  "09111 Xenobiotics biodegradation and metabolism",
  "09191 Unclassified: metabolism",
  "09194 Poorly characterized",
  "09193 Unclassified: signaling and cellular processes",
  "09192 Unclassified: genetic information processing",
  "09145 Cellular community - prokaryotes"
)

clean_label    = function(x) gsub("^\\d+ ", "", gsub(" \\[.*?\\]", "", x))
binom_p        = function(k, n, p, alt = "greater") binom.test(k, n, p = p, alternative = alt)$p.value
save_plot      = function(path, p, ...) ggsave(file.path(OUT, path), p, ...)
load_rds       = function(name) readRDS(file.path(DAT, paste0(name, ".rds")))

# duplicate level c names 
disambiguate_C = function(dat) {
  dat %>% mutate(label_C = case_when(
    Level_C == "03050 Proteasome [PATH:ko03050]"  ~ "Proteasome (pathway)",
    Level_C == "03051 Proteasome [BR:ko03051]"    ~ "Proteasome (BRITE)",
    Level_C == "03010 Ribosome [PATH:ko03010]"    ~ "Ribosome (pathway)",
    Level_C == "03011 Ribosome [BR:ko03011]"      ~ "Ribosome (BRITE)",
    Level_C == "03040 Spliceosome [PATH:ko03040]" ~ "Spliceosome (pathway)",
    Level_C == "03041 Spliceosome [BR:ko03041]"   ~ "Spliceosome (BRITE)",
    TRUE ~ label_C
  ))
}
