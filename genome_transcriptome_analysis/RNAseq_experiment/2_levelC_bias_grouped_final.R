# level c pathway up/downregulation

# output dir
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# data
source("genome_transcriptome_analysis/RNAseq_experiment/utils.R")
sig_by_pathway = load_rds("sig_by_pathway")
null_props     = load_rds("null_props")
p_up_global    = null_props$p_up
p_down_global  = null_props$p_down

# per pathway stats
pathway_bias = sig_by_pathway %>%
  group_by(Level_B, Level_C) %>%
  summarise(
    n_up     = sum(log2FoldChange > 0),
    n_down   = sum(log2FoldChange < 0),
    n_total  = n(),
    mean_lfc = mean(log2FoldChange),
    .groups  = "drop"
  ) %>%
  filter(n_total >= 10) %>%
  mutate(
    pct_up   = n_up   / n_total * 100,
    pct_down = n_down / n_total * 100,
    net_bias = pct_up - pct_down,
    label_C  = clean_label(Level_C),
    label_B  = clean_label(Level_B)
  ) %>%
  disambiguate_C()

# binomial tests w/ bh correction
pathway_bias = pathway_bias %>%
  rowwise() %>%
  mutate(
    p_down_raw = binom_p(n_down, n_total, p_down_global),
    p_up_raw   = binom_p(n_up,   n_total, p_up_global)
  ) %>%
  ungroup()

n_paths  = nrow(pathway_bias)
padj_all = p.adjust(c(pathway_bias$p_down_raw, pathway_bias$p_up_raw), method = "BH")

pathway_bias = pathway_bias %>%
  mutate(
    padj_down = padj_all[seq_len(n_paths)],
    padj_up   = padj_all[seq_len(n_paths) + n_paths],
    sig_down  = padj_down < 0.05,
    sig_up    = padj_up   < 0.05,
    sig_any   = sig_up | sig_down
  )

# output 
pathway_bias %>%
  select(label_B, label_C, n_total, n_up, n_down, pct_up, pct_down, net_bias,
         mean_lfc, padj_down, padj_up, sig_down, sig_up) %>%
  arrange(desc(net_bias)) %>%
  write.csv(file.path(OUT, "levelC_bias_stats.csv"), row.names = FALSE)

# pathways to plot
# all significantly pathways (padj < 0.05)
# trending pathways: >75% in one direction AND n >= 20
plot_paths = pathway_bias %>%
  filter(
    sig_any |                                    
    (pct_up > 75 & n_total >= 20) |          
    (pct_down > 75 & n_total >= 20)           
  ) %>%
  arrange(desc(mean_lfc))

# assign func categories 
assign_category = function(levelC_name, levelB_name) {
  case_when(
    str_detect(levelC_name, regex("chaperone|protein folding|ubiquitin|proteasome|protein processing|unfolded protein",
                                  ignore_case = TRUE)) ~ "Proteostasis",
    str_detect(levelC_name, regex("DNA repair|homologous recombination|mismatch repair|base excision|nucleotide excision|replication|fanconi anemia|non-homologous end joining|NHEJ|double strand break|chromosome|condensin|cohesin|chromatin|nucleosome|histone",
                                  ignore_case = TRUE)) ~ "DNA Repair & Genome Maintenance",
    str_detect(levelC_name, regex("biosynthesis|lipid biosynthesis|starch|sucrose|sphingolipid|peptidase",
                                  ignore_case = TRUE)) ~ "Biosynthesis",
    str_detect(levelC_name, regex("glycolysis|TCA|electron transport|oxidative phosphorylation|glyoxylate|dicarboxylate|amino acid|carbon|calvin|pentose phosphate|glutathione|purine metabolism|folate|valine|leucine|isoleucine|cysteine|methionine|glycine|serine|threonine|alanine|aspartate|glutamate",
                                  ignore_case = TRUE)) ~ "Energy Metabolism",
    str_detect(levelC_name, regex("RNA processing|splicing|transcription|ribosome biogenesis|RNA helicases|nucleocytoplasmic transport",
                                  ignore_case = TRUE)) ~ "RNA Processing & Gene Expression",
    str_detect(levelC_name, regex("SNARE|vesicle|trafficking|endoplasmic reticulum|golgi|COPI|COPII|autophagy",
                                  ignore_case = TRUE)) ~ "Membrane Trafficking & Vesicles",
    str_detect(levelC_name, regex("cell cycle|mitosis|meiosis|cytokinesis",
                                  ignore_case = TRUE)) ~ "Cell Cycle & Division",
    str_detect(levelC_name, regex("actin|tubulin|cytoskeleton|kinesin|dynein|motor protein",
                                  ignore_case = TRUE)) ~ "Cytoskeleton & Motility",
    str_detect(levelC_name, regex("signaling|MAPK|ras|calcium|receptor|hormone",
                                  ignore_case = TRUE)) ~ "Cell Signaling Pathways",
    str_detect(levelC_name, regex("kinase|phosphatase|phosphorylation",
                                  ignore_case = TRUE)) ~ "Protein Kinases & Phosphorylation",
    
    # others from lvl b
    str_detect(levelB_name, regex("Genetic Information Processing", ignore_case = TRUE)) ~ "RNA Processing & Gene Expression",
    str_detect(levelB_name, regex("Environmental Information Processing", ignore_case = TRUE)) ~ "Cell Signaling Pathways",
    str_detect(levelB_name, regex("Metabolism", ignore_case = TRUE)) ~ "Energy Metabolism",
    
    TRUE ~ "Other"
  )
}

plot_paths = plot_paths %>%
  mutate(
    category = assign_category(label_C, label_B),
    category = factor(category, levels = c(
      "Proteostasis", 
      "DNA Repair & Genome Maintenance", 
      "RNA Processing & Gene Expression", 
      "Membrane Trafficking & Vesicles", 
      "Cell Cycle & Division",
      "Cytoskeleton & Motility", 
      "Cell Signaling Pathways", 
      "Protein Kinases & Phosphorylation",
      "Energy Metabolism", 
      "Biosynthesis", 
      "Other"
    ))
  ) %>%
  arrange(category, desc(mean_lfc)) %>%
  mutate(label_C = factor(label_C, levels = rev(label_C)))  # reverse for ggplot y-axis

# plot prep
plot_long = plot_paths %>%
  transmute(
    label_C, label_B, category, n_total, sig_down, sig_up,
    pct_down = -pct_down,     # negative for leftward bars
    pct_up
  ) %>%
  pivot_longer(c(pct_down, pct_up), names_to = "direction", values_to = "pct") %>%
  mutate(
    direction = recode(direction, pct_down = "Down", pct_up = "Up"),
    is_sig    = case_when(
      direction == "Down" & sig_down ~ TRUE,
      direction == "Up"   & sig_up   ~ TRUE,
      TRUE ~ FALSE
    ),
    # sig pathways vs trending
    alpha_val = if_else(is_sig, 0.88, 0.35),
    label_C   = factor(label_C, levels = levels(plot_paths$label_C))
  )

# labels
n_labels = plot_paths %>%
  mutate(
    x_pos = pct_up + 1.5,
    label = paste0("(n=", n_total, ")")
  )

# breaks
category_breaks = plot_paths %>%
  group_by(category) %>%
  summarise(
    first_pos = which(levels(plot_paths$label_C) == last(label_C)),  # last because y-axis is reversed
    last_pos = which(levels(plot_paths$label_C) == first(label_C)),
    mid_pos = (first_pos + last_pos) / 2,
    .groups = "drop"
  ) %>%
  filter(category != "Other")  # Don't show "Other" category headers


# plot
p_comprehensive = ggplot(plot_long,
    aes(y = label_C, x = pct, fill = direction, alpha = alpha_val)) +

  geom_col(position = "identity", width = 0.75) +

  # n= count labels
  geom_text(data = n_labels,
            aes(y = label_C, x = x_pos, label = label),
            size = 2.8, colour = "grey50", hjust = 0, inherit.aes = FALSE) +

  # Category headers on the left (outside plot area)
  geom_text(data = category_breaks,
            aes(y = mid_pos, x = -125, label = category),
            size = 3.2, colour = "grey30", fontface = "bold", 
            hjust = 0, inherit.aes = FALSE) +

  # Category separators
  geom_hline(data = category_breaks %>% slice(-n()),  # exclude last separator
             aes(yintercept = first_pos - 0.5), 
             colour = "grey80", linewidth = 0.3, inherit.aes = FALSE) +

  # zero reference line
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60") +

  scale_fill_manual(
    values = c(Down = COL_DOWN, Up = COL_UP),
    labels = c(Down = "% downregulated", Up = "% upregulated")
  ) +
  scale_alpha_identity() +

  scale_x_continuous(
    limits = c(-130, 115),
    breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100),
    labels = function(x) paste0(abs(x), "%"),
    expand = c(0, 0)
  ) +

  labs(
    x = expression(symbol("\254") ~ "% downregulated" ~ ~ "% upregulated" ~ symbol("\256")),
    y = NULL,
    fill = NULL
  ) +

  theme_rna() +
  theme(
    legend.position = "top",
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.text.y     = element_text(size = 7),
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 10, colour = "grey40"),
    plot.margin     = margin(8, 12, 8, 45) 
  )

p_comprehensive
