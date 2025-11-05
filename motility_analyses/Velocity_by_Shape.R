#bin velocity by aspect ratio category this is for one temp only 
library(dplyr)

Verm_df <- data[data$temp == "25V", ]
high_df<-  data[data$temp == "57", ]
df <- Verm_df %>%
  mutate(
    shape = case_when(
      #AR >= 1 & AR < 1.25 ~ "cyst",
      AR >= 1.0 & AR < 2.5 ~ "lamelli",
      AR >= 2.5 & AR < 20 ~ "worm",
      FALSE ~ "Other"
    )
  )

#find runs of X frames where shape doesn't change and capture ave velocity 
df_runs <- df %>%
  arrange(group, ID, Frame) %>%
  group_by(group, ID) %>%
  mutate(run_id = with(rle(shape), rep(seq_along(lengths), lengths))) %>%
  group_by(group, ID, run_id) %>%
  filter(n() >= 3) %>%
  summarise(
    shape = first(shape),
    avg_speed = mean(stepwise_speed * 60, na.rm = TRUE),
    avg_AR = mean(AR, na.rm = TRUE),  # <-- add this line
    frames_in_run = n(),
    temp = first(temp),
    .groups = "drop"
  )
#View(df_runs)
speed_by_shape <- df_runs %>%
  group_by(shape) %>%
  summarise(
    median_velocity = median(avg_speed, na.rm = TRUE),
    sd_velocity = sd(avg_speed, na.rm = TRUE),
    n_runs = n(),
    sem_velocity    = sd_velocity / sqrt(n_runs),
    .groups = "drop"
  )
#View(speed_by_shape)
write.csv(speed_by_shape, "shape_velocity_25V.csv")


switching_rate <- df_runs %>%
  arrange(group, ID, run_id) %>%
  group_by(group, ID) %>%
  mutate(shape_change = shape != lag(shape)) %>%
  summarise(
    n_switches = sum(shape_change, na.rm = TRUE),
    total_time = sum(frames_in_run) / 60,  # convert to minutes if frames are per second
    rate_per_min = n_switches / total_time,
    #sem_switch    = sd_velocity / sqrt(n_runs),
    temp = first(temp),
    .groups = "drop"
  )
View(switching_rate)

shape_speed_summary <- switching_rate %>%
  group_by(temp) %>%
  summarise(
    switchrate = mean(rate_per_min, na.rm = TRUE),
    #runlength = mean(ave_run_length, na.rm = TRUE),
    #sd_speed = sd(avg_speed, na.rm = TRUE),
    #n_segments = n(),
    #se = sd_speed / sqrt(n_segments),
    #ci_low = speed - 1.96 * se,
    #ci_high = speed + 1.96 * se,
    #.groups = "drop"
  )
View(shape_speed_summary)
#write.csv(shape_speed_summary, "shape_velocity_binned_25V.csv")

switching_summary <- df_runs %>%
  group_by(group,ID) %>%
  summarise(
    n_runs = n(),
    total_time = sum(frames_in_run),
    switches = n_runs - 1,
    switch_rate_per_min = (switches*60 / (total_time / 1)), # set frame_rate if not 1
    avg_run_length = mean(frames_in_run),
    .groups = "drop"
  )
View(switching_summary)

p <- ggplot(df_runs, aes(x = shape, y = avg_speed, group = shape)) +
  coord_cartesian(ylim = c(0, 50)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(
    dodge.width = 1,
    aes(color = shape),
    size = 0.5,
    alpha = 0.5,
    cex = 0.25
  ) +
  scale_color_manual(values = custom_colors)

p <- p + theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1) # keeps shape labels readable
  ) +
  ylab("Average Speed (µm/min)") +
  xlab("Shape Category")

p

write.csv(switching_summary, "shape_switch_binned_25V.csv")


p<-ggplot(data, aes(x = factor(temp), y = speed, color = shape, group = shape)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  coord_cartesian(ylim = c(10, 35)) +
  #geom_errorbar(aes(ymin = speed - sd_speed, ymax = speed + sd_speed),
              #  position = position_dodge(width = 0.3), width = 0.2) +
  labs(x = "Temperature (°C)", y = "Speed (µm/min)", color = "Shape") +
  theme_classic(base_size = 14)
p

# switching ---------------------------------------------------------------
#n_runs: How many continuous shape segments exist for that ID.
#avg_run_length_sec: Average time (sec) spent in one shape before switching.
#amelli / worm / Other: Total time (sec) spent in each shape.
library(dplyr)
library(tidyr)

# ---- Shape classification ----
df <- Verm_df %>%
  mutate(shape = case_when(
    AR >= 1.1 & AR < 2.5 ~ "lamelli",
    AR >= 2.5 & AR < 20  ~ "worm"
  )) %>%
  filter(!is.na(shape))  # is shape NA?

frame_rate <- 1  # make sure this frame rate is right####################

# run IDs for continuous shape
df_runs <- df %>%
  arrange(group, ID, Frame) %>%
  group_by(group, ID) %>%
  mutate(run_id = cumsum(shape != lag(shape, default = first(shape)))) %>%
  group_by(group, ID, run_id) %>%
  filter(n() >= 5) %>%  # keep only runs of at least 5 frames, adjust this number
  summarise(
    shape = first(shape),
    avg_speed = mean(stepwise_speed * 60, na.rm = TRUE),  # um/min if speed in um/frame otherwise change this
    frames_in_run = n(),
    temp = first(temp),
    .groups = "drop"
  )

# summarize for metrics 
switching_summary <- df_runs %>%
  group_by(group, ID) %>%
  summarise(
    n_runs = n(),
    switches = n_runs - 1,  ##### switches = runs - 1 ########
    total_frames = sum(frames_in_run),
    switch_rate_per_min = (switches * 60 * frame_rate) / total_frames,  # switches per minute
    avg_run_length_sec = mean(frames_in_run) / frame_rate,
    temp = first(temp),
    .groups = "drop"
  )

# ---- Total time spent in each shape (wide format) ----
shape_time <- df_runs %>%
  group_by(group, ID, shape) %>%
  summarise(time_in_shape_sec = sum(frames_in_run) / frame_rate, .groups = "drop") %>%
  pivot_wider(names_from = shape, values_from = time_in_shape_sec, values_fill = 0)

# combine it
final_summary <- switching_summary %>%
  left_join(shape_time, by = c("group", "ID"))

print(final_summary)
name <- "25V" 
#final_summary$temp <- 60
write.csv(df, paste0(name, "_shapefreq.csv"))



# Calculate average switch_rate_per_min per group
ave_switch_rate <- final_summary %>%
  mutate(interval_per_min = 1 / switch_rate_per_min) %>%  # create the new column
  #group_by(group) %>%
  summarise(
    median_switch_rate_per_min = median(switch_rate_per_min, na.rm = TRUE),
    sd_switch_rate_per_min = sd(switch_rate_per_min, na.rm = TRUE),
    n = sum(!is.na(switch_rate_per_min)),
    sem_switch_rate_per_min = sd_switch_rate_per_min / sqrt(n),
    mean_interval_per_min = mean(interval_per_min, na.rm = TRUE),
    .groups = "drop"
  )

print(ave_switch_rate)
write.csv(ave_switch_rate,"aveshapefreq_Verm-sem.csv")


#get counts
verm_unique <- Verm_df %>%
  group_by(group) %>%
  summarise(unique_particles = n_distinct(ID)) %>%
  summarise(total = sum(unique_particles)) %>%
  pull(total)



