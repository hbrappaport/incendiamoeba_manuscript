# ordination of 48 vs. 61C 
source("genome_transcriptome_analysis/RNAseq_experiment/utils.R")

# data
all_genes = load_rds("all_genes")

# all genes
mat_all = all_genes %>%
  select(all_of(SAMPLE_COLS)) %>%
  as.matrix()

mat_all = mat_all - rowMeans(mat_all)
pca_all  = prcomp(t(mat_all), scale. = FALSE)
pct_all  = round(pca_all$sdev^2 / sum(pca_all$sdev^2) * 100, 1)

scores_all = as_tibble(pca_all$x[, 1:2]) %>%
  mutate(sample = SAMPLE_COLS) %>%
  left_join(SAMPLE_INFO, by = "sample")

p_pca_all = ggplot(scores_all, aes(x = PC1, y = PC2,
                                   fill = condition,
                                   label = sample)) +
  geom_point(size = 5, shape = 21, colour = "black",
             stroke = 0.5, alpha = 0.95) +
  scale_fill_manual(values = c("48C" = COL_DOWN,
                               "61C" = COL_UP),
                    name = "Condition") +
  labs(x = paste0("PC1 (", pct_all[1], "%)"),
       y = paste0("PC2 (", pct_all[2], "%)")) +
  theme_classic() + coord_equal()

# outputs
p_pca_all
