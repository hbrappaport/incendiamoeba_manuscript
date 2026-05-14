# packages
library(httr)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

# TPM gene level from RNA seq at 48ºC and 61ºC
tpms_gene = read.csv("genome_transcriptome_analysis/RNAseq_experiment/results_pilot3_w_KO_PFAM_TPMs.csv") 
tpms_gene = tpms_gene %>%
  distinct(Geneid, .keep_all = TRUE) # remove dupes from kegg annot

# select top 500 most abundant genes at 61ºC (by transcripts per million)
top500_61gene = tpms_gene %>%
  arrange(desc(average_TPM_61)) %>%
  slice(1:500)
top500_genelist = top500_61gene %>%
  select(Geneid)
write.csv(top500_genelist, "genome_transcriptome_analysis/RNAseq_experiment/top500_genelist_61average.csv")

# % of transcriptional activity of the top 500, rel abundance of top 500 averaged across samples
# total transcripts across samples 61ºC
sum_transcripts = sum(tpms_gene$average_TPM_61)
# add column for relative abundance of each gene
tpms_relabund = tpms_gene %>%
  mutate(relabund = average_TPM_61/sum_transcripts)
# sum relative abundance of top 500 genes
top500_61relabund = tpms_relabund %>%
  arrange(desc(relabund)) %>%
  slice(1:500)
sum_top500 = sum(top500_61relabund$relabund) # 49.9%

# fraction of transcripts by gene + cutoff at 500
tpms_accumulate = tpms_relabund %>%
  mutate(accumulate = cumsum(relabund)) 
ggplot(tpms_accumulate, aes(reorder(Geneid, accumulate), accumulate, fill = "Geneid")) +
  geom_bar(stat="identity") +
  #scale_y_log10() +
  scale_fill_manual(values = c("#C1532D")) +
  geom_vline(xintercept = 500, linetype = "dotted", size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", size = 1) +
  theme_classic() +  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position="none") +
  labs(y ="% of Transcripts", x= "Gene")

# select most upregulated
top50up_61gene = tpms_gene %>%
  arrange(desc(log2FoldChange)) %>%
  slice(1:50)
top50up_genelist = top50up_61gene %>%
  select(Geneid)
write.csv(top50up_genelist, "genome_transcriptome_analysis/RNAseq_experiment/top50up_genelist.csv")

# select most downregulated
top50down_61gene = tpms_gene %>%
  arrange(log2FoldChange) %>%
  slice(1:50)
top50down_genelist = top50down_61gene %>%
  select(Geneid)
write.csv(top50down_genelist, "genome_transcriptome_analysis/RNAseq_experiment/top50down_genelist.csv")
