# packages
library(tidyverse)

# data
tms = read.csv("protein_analyses/4_Tm_analysis/allTms.csv")

# filter to just Vermamoeba & one Incendiamoeba
tms_filt = tms %>%
  filter(species_braker_file %in% c("Vermamoeba_vermiformisCDC19",
                                    "Incendiamoeba_cascadensis2"))

ggplot(tms_filt, aes(Tm, fill=species_braker_file)) +
  geom_density(alpha=0.75) +
  scale_fill_manual(values = c("#C1532D",
                               "#d1d1d1")) +
  geom_vline(xintercept = 50.9, linetype = "dashed") +
  geom_vline(xintercept = 46.9, linetype = "dashed") +
  theme_classic()


tms_sum = tms %>%
  group_by(species_braker_file) %>%
  summarise(median = median(Tm),
            mean = mean(Tm),
            total = sum(Tm),
            above_60_pct = sum(Tm>60)/sum(Tm)*100)

# filter tms above 60ºC
tms_above_60 = tms %>%
  filter(Tm >= 60)

# plot tms above 60ºC
ggplot(tms_above_60, aes(species_braker_file, Tm)) +
  geom_jitter() +
  geom_boxplot()

