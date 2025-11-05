# packages
library(ggplot2)
library(dplyr)
library(gcplyr)
library(tidyr)
library(broom)

# data
gcurves = read.csv("growth_doubling/growthcurves_individual.csv")
doublingtime = read.csv("growth_doubling/incendi_doublingtimes.csv")

# Growth curves

# color pal
temp_colors = c(
  "30C"   = "#313695",
  "40C"   = "#4575B4",
  "42C"   = "#9FCAFF",
  "45C"   = "#A6D96A",
  "50C"   = "#FFFFD4",
  "51.5C" = "#FFE600",
  "54C"   = "#FFC300",
  "55C"   = "#D3BB9C",
  "57C"   = "#FDAE61",
  "61C"   = "#F46D43",
  "61.5C" = "#E34A33",
  "62C"   = "#D73027",
  "62.5C" = "#A50026",
  "63C"   = "#CD782E",
  "64C"   = "#C1532D"
)

# plot growth curves at each incubation temp
ggplot(gcurves, aes(Time, Cell_mL)) +
  geom_smooth(aes(color=TempC), span=3,
              se = T, linewidth=0.75,
              alpha = 0.1) +
  geom_jitter(aes(fill=TempC), pch=21, stroke=0.2, width = 0.01) +
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values = temp_colors) +
  ylab("Cells/mL") +
  theme_classic(base_size = 14) +
  facet_wrap(~ TempC, scales = "free")


# log transform
gcurves$log_cells = log(gcurves$Cell_mL)

# calc summary stats
growth_summary = gcurves %>%
  group_by(Batch_ID) %>%
  do(tidy(lm(log_cells ~ Time, data = .))) %>%
  filter(term == "Time") %>%
  mutate(sig = ifelse(p.value < 0.05 & estimate > 0,
                      "Yes", "No"))

# Doubling time

# filter to most representative batch, considering temp probe measurements
doublingtime_best = doublingtime %>% filter(doublingtime$Best_batch == TRUE)

# plot doubling time curve
ggplot(doublingtime_best, aes(x = Temperature, y = Average.Doubling.Time..Exponential.portion_hours.)) +
  geom_smooth(span = 1, fill= "#D1D1D1", color = "#000000") +
  geom_point(size = 1) +
  scale_y_continuous(breaks = seq(25, 125, by = 25)) +
  ylab("Average doubling time (hours)") +
  xlab("Temperature (ºC)") +
  theme_classic(base_size = 14)

