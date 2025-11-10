# libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

genomestats = read.csv("genome_transcriptome_analysis/genomestats.csv")
busco = read.csv("genome_transcriptome_analysis/busco.csv")
orthocompare = read.csv("genome_transcriptome_analysis/orthogroup_comparison.csv")

# Busco Completeness/Contamination
ggplot(busco, aes(fill=Complete.Contam, y=Value, x=Genome)) + 
  geom_bar(position= position_stack(reverse=TRUE), stat="identity") +
  scale_fill_manual(values = c("#5389c9", "#e8881a")) +
  coord_flip() +
  theme_classic()


# Total length
ggplot(genomestats, aes(fill = "#9fcaff" , y = Total.length, x = Genome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#9fcaff")) +
  coord_flip() +
  theme_classic()

# Percent GC
ggplot(genomestats, aes(fill = "#9fcaff", y = X..GC, x = Genome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#c3ddff")) +
  coord_flip() +
  theme_classic()

# Contigs
ggplot(genomestats, aes(fill = "#9fcaff", y = X..Contigs, x = Genome)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#5389c9")) +
  coord_flip() +
  theme_classic()

#Orthogroup Comparison
orthocompare_longer = pivot_longer(orthocompare, !species)
# tiles for dram functions
ggplot(orthocompare_longer) +
  geom_tile( 
    aes(x = species, y = name, fill = value)) +
  scale_fill_gradient(
    low = "white",
    high = "#C1532D"
  ) +
  # scale_color_manual(values = "white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  coord_flip()
