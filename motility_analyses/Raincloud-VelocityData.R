#raincloud plots, velocity data
library(ggplot2)
library(ggdist)       # cloud
library(dplyr)
library(ggpp)
custom_colors <- c(
  "#5662B2",
  "#6389c9",
  "#9FCAFF",
  "#d1d1d1",
  "#D3BB9C",
  "#D7A465",
  "#F09900",
  "#e9881a",
  "#CD782E",
  "#C1532D",
  "#AD372A"
)
gray <- c(
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0",
  "#C0C0C0"
)
# raincloud plot mean velocity per group or per group,ID ------------------
# data structure: #   temperature #   IDs individual track #group - movie
# stepwise_speed — frame-by-frame speed; averaged per group
# 25V=vermamoeba vermiformis
p <- ggplot(mapping = aes(x = temp, y = velocity,
                          fill = temp, color = temp)) +
#id level
  stat_halfeye(
    data         = ID_velocity %>% filter(temp != "25V"),
    adjust       = 0.5,
    width        = 0.6,
    .width       = 0,
    point_colour = NA,
    slab_alpha   = 0.6,
    position     = position_nudge(x = 0.15)
  ) +
  geom_boxplot(
    data          = group_velocity %>% filter(temp != "25V"),
   width         = 0.4,
    outlier.shape = NA,
   color         = "black",
   fill          = NA,
   position      = position_nudge(x = 0.35)
 ) +
  # group level
  stat_dots(
    data      = group_velocity %>% filter(temp != "25V"),
    side      = "left",
    dotsize   = 0.9,
    binwidth  = 0.01,
    alpha     = 7,
    position  = position_nudge(x = -0.05)
  ) +
  scale_fill_manual(values  = gray) +
  scale_color_manual(values = gray) +
  coord_cartesian(ylim = c(0, 40)) +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  labs(
    y     = "Velocity (um/min)",
    x     = "Temperature"
  ) +
  theme_classic() +
  theme(legend.position = "none")
p


