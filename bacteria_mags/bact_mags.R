# packages
library(ggplot2)
library(dplyr)
library(tidyr)

magstats = read.csv("bacteria_mags/bestmags_bowtie.csv")
dram_lassen = read.csv("bacteria_mags/lassen_metag_withetc.csv")
dram_lassen_longer = pivot_longer(dram_lassen, cols = !genome)

# bars for MAG abundance
ggplot(magstats, aes(y = GTDB.Top.hit.Species, x = Bowtie.mapped.reads...from.merged.read.set)) +
  geom_col(fill = "#5388c9") +
  geom_text(
    aes(label = Bowtie.mapped.reads...from.merged.read.set), 
    hjust = -0.1
  ) +
  theme_classic()

# bars for complete/contam
ggplot(magstats, aes(y = GTDB.Top.hit.Species, x = X..Complete)) +
  geom_col(fill = "#5388c9") +
  theme_classic()

ggplot(magstats, aes(y = GTDB.Top.hit.Species, x = X..Contam)) +
  geom_col(fill = "#5388c9") +
  theme_classic()

# tiles for dram functions
ggplot(dram_lassen_longer) +
  geom_tile( 
    aes(x = name, y = genome, fill = value)) +
  scale_fill_gradient(
    low = "lightgrey",
    high = "#bd6e15"
  ) +
  # scale_color_manual(values = "white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

