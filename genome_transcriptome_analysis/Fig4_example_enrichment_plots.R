## Examples of IND KOs enriched 
library(tidyverse)

# data
kb = read.csv("genome_transcriptome_analysis/keggbrite_5genomes.csv")

# look at major categories
lc = as.data.frame(sort(table(kb$Level_C)))
table(kb$Level_B)

kb_cats = kb %>%
  select(Level_C, Level_B) %>%
  distinct()

# filter to remove non focal categories
kb = kb %>%
  filter(!Level_B %in% c("09176 Drug resistance: antineoplastic",
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
                         "09145 Cellular community - prokaryotes"))

#################################################################
#################################################################

# plots looking at individual KOs
kb_long = kb %>%
  select(-Level_B, - Level_C) %>%
  distinct() %>%
  select(KO_ID, Description, Avg_excess,
         vverm, acastellanii, dicty,
         bc2001, bc2002) %>%
  pivot_longer(
    cols = -c(KO_ID, Description, Avg_excess),
    names_to = "genome",
    values_to = "gene_copy_detected"
  )


# HSPs
HSP20 = ggplot(kb_long %>% filter(KO_ID == "K13993"),
               aes(fct_reorder(genome,-gene_copy_detected),
                   gene_copy_detected)) +
  geom_bar(stat="identity", aes(fill= genome),
           color="black", size = .5) +
  scale_fill_manual(values=c("gray", "#C1532D", "#C1532D", "gray", "gray")) +
  theme_classic() +
  labs(x = "Genome", y = "Gene Copy",
       title = "K13993") +
  theme(legend.position="none")

# Ubiquitin-protein ligases - RLIM
RLIM = ggplot(kb_long %>% filter(KO_ID == "K16271"),
              aes(fct_reorder(genome,-gene_copy_detected),
                  gene_copy_detected)) +
  geom_bar(stat="identity", aes(fill= genome),
           color="black", size = .5) +
  scale_fill_manual(values=c("gray", "#C1532D", "#C1532D", "gray", "gray")) +
  theme_classic() +
  labs(x = "Genome", y = "Gene Copy",
       title = "K16271") + 
  theme(legend.position="none")

# Zinc fingers -  ZBED1; zinc finger BED domain-containing protein 1 (E3 SUMO-protein ligase ZBED1) [EC:2.3.2.-]
ZBED1 = ggplot(kb_long %>% filter(KO_ID == "K24637"),
               aes(fct_reorder(genome,-gene_copy_detected),
                   gene_copy_detected)) +
  geom_bar(stat="identity", aes(fill= genome),
           color="black", size = .5) +
  scale_fill_manual(values=c("gray", "#C1532D", "#C1532D", "gray", "gray")) +
  theme_classic() +
  labs(x = "Genome", y = "Gene Copy",
       title = "K24637") + 
  theme(legend.position="none")

### plot all together
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## And plot
multiplot(ZBED1, RLIM, HSP20,cols = 3)



