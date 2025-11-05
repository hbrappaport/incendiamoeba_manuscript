#mean squared displacement
library(dplyr)
library(broom)
library(tidyr)
library(purrr)

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


#you want to combine all the output files for analysis into one file: 
#first set WD to a folder with all the individual files to combine 
path <- getwd()
csvcomb <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
datacomb <- bind_rows(lapply(csvcomb, read_csv))
# Save the combined data to a new CSV file, rename it 
write_csv(datacomb, "filename.csv")
#open it and check that the files are now combined so that they all share column names
#reimport the summary combined data as df 
df<-data
#use your dataframe of raw combined outputs

# How many groups per temp before filtering?
df %>%
  group_by(temp, group, ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(temp, name = "n_tracks_before")

# recode ------------------------------------------------------------------
# filter and map 
df<-data
df <- df %>% filter(stepwise_displacement < 10) #anything over this is detached

# How many groups per temp before filtering?
df %>%
  group_by(temp, group, ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(temp, name = "n_tracks_before")
#is this not working? check what group is assigned as (it should be each individual movie)

# Compute the dynamic max_tau as 1/4th of the shortest track length
max_tau <- floor(median(df %>% group_by(group,ID) %>% summarise(n = n()) %>% pull(n)))
#max_tau <- floor(min(df %>% group_by(group, ID) %>% summarise(n = n()) %>% pull(n)))
#max_tau : the maximum time lag over which you calculate displacements

valid_groups <- df %>%
  group_by(group, temp,ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > max_tau)

msd_df <- df %>%
  group_by(group, temp, ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > max_tau) %>%
  mutate(group_id = paste(group, ID, sep = "_")) %>%
  pmap_dfr(function(group, temp, ID, n, group_id) {
    track <- df %>%
      filter(group == !!group, ID == !!ID) %>%
      arrange(Frame)
  #max_tau is computed based on track length to avoid bias in particle tracking 
    tibble(
      group = group,
      temp = temp,
      ID = ID,
      tau = 1:max_tau,
      MSD = map_dbl(1:max_tau, function(tau) {
        disp <- track %>%
          mutate(
            X_shift = lead(Xpos, tau),
            Y_shift = lead(Ypos, tau)
          ) %>%
          filter(!is.na(X_shift)) %>%
          mutate(sq_disp = (X_shift - Xpos)^2 + (Y_shift - Ypos)^2) %>%
          summarise(msd = mean(sq_disp, na.rm = TRUE)) %>%
          pull(msd)
        return(disp)
      })
    )
  })

msd_overall <- msd_df %>%
  group_by(temp, tau) %>%
  summarise(
    MSD_mean = mean(MSD, na.rm = TRUE),
    MSD_se = sd(MSD, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(time = tau / max(tau))  # Normalize to 0–1, only works if time interval is 1s which it is 

df %>%
  group_by(temp, group,ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  count(temp, name = "tracks_per_temp")

# net displacement --------------------------------------------------------
df<-data
# Compute net displacement for each track (group + ID)
net_disp_df <- df %>%
  group_by(temp, group, ID) %>%
  summarise(
    start_x = first(Xpos),
    start_y = first(Ypos),
    end_x = last(Xpos),
    end_y = last(Ypos),
    net_displacement = sqrt((end_x - start_x)^2 + (end_y - start_y)^2),
    .groups = "drop"
  )

# by temp
net_disp_summary <- net_disp_df %>%
  group_by(temp) %>%
  summarise(
    mean_net_disp = mean(net_displacement, na.rm = TRUE),
    median_net_disp = median(net_displacement, na.rm = TRUE),
    sd_net_disp = sd(net_displacement, na.rm = TRUE),
    n_tracks = n(),
    .groups = "drop"
  )

# summary table
net_disp<-as.data.frame(net_disp_summary)
write.csv(net_disp, "net_displacement.csv")

# straightness ratio: net displacement, path length, and straightness per track
disp_metrics_df <- df %>%
  group_by(temp, group, ID) %>%
  summarise(
    start_x = first(Xpos),
    start_y = first(Ypos),
    end_x = last(Xpos),
    end_y = last(Ypos),
    net_displacement = sqrt((end_x - start_x)^2 + (end_y - start_y)^2),
    path_length = sum(sqrt(diff(Xpos)^2 + diff(Ypos)^2), na.rm = TRUE),
    straightness = ifelse(path_length > 0, net_displacement / path_length, NA),
    .groups = "drop"
  )

# Compute reference mean straightness 
reference_mean <- disp_metrics_df %>%
  filter(temp %in% c(25, 70)) %>%
  summarise(ref_mean = mean(straightness, na.rm = TRUE)) %>%
  pull(ref_mean)

# Normalize straightness
disp_metrics_df <- disp_metrics_df %>%
  mutate(norm_straightness = straightness / reference_mean)

# normalize summ
straightness_summary <- disp_metrics_df %>%
  group_by(temp) %>%
  summarise(
    mean_straightness = mean(straightness, na.rm = TRUE),
    median_straightness = median(straightness, na.rm = TRUE),
    mean_norm_straightness = mean(norm_straightness, na.rm = TRUE),
    median_norm_straightness = median(norm_straightness, na.rm = TRUE),
    n_tracks = n(),
    .groups = "drop"
  )
View(straightness_summary)
#out
write.csv(disp_metrics_df, "raw_normalized_straightness.csv", row.names = FALSE)
write.csv(straightness_summary, "straightness.csv", row.names = FALSE)

# plot --------------------------------------------------------------------
p<-ggplot(msd_overall, aes(x = time, y = MSD_mean, color = as.factor(temp))) +
  geom_line(size = 1) + 
  ylim(0,750)+
  geom_ribbon(aes(ymin = MSD_mean - MSD_se, ymax = MSD_mean + MSD_se, fill = as.factor(temp)), alpha = 0.08) +
  geom_point(size = 0) +
  labs(
    title = "MSD Comparison Across Temperatures",
    x = "Time (seconds)",
    y = "Mean Squared Displacement (MSD)",
    color = "Temperature",
    fill = "Temperature"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  #scale_color_viridis_d() +  # Viridis color scale for lines
  #scale_fill_viridis_d() +   # Viridis color scale for ribbons
  theme_minimal() +                          
  theme(
    text = element_text(size = 14),
    panel.grid = element_blank(),    # Remove gridlines
    axis.line = element_line(color = "black"),  # Keep axis lines
    legend.title = element_blank()   # Remove legend title
  ) 
p


#plotting via facet wrap, this is for multiple groups but you should only have one group 
#change the code so that it shows you just the one group that you ran those videos for
p <- ggplot(msd_overall, aes(x = time, y = MSD_mean, color = as.factor(temp))) +
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = MSD_mean - MSD_se, ymax = MSD_mean + MSD_se, fill = as.factor(temp)), alpha = 0.2) +
  geom_point(size = 0) +
  labs(
    title = "MSD Comparison Across Temperatures",
    x = "Time (seconds)",
    y = "Mean Squared Displacement (MSD)",
    color = "Temperature",
    fill = "Temperature"
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  #scale_color_viridis_d() +
  #scale_fill_viridis_d() +
  facet_wrap(~ temp, ncol = 4, scales = "free") +
  coord_cartesian(ylim = c(0, 600)) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.title = element_blank()
  )

p

#save this plot, how does it compare to the plots in the manuscript? 

#just read through some of these parts to get a general idea of what's happening
# Analysis::  -------------------------------------------------------------
#get slopes log-transform tau and MSD, filter invalid points
msd_slopes <- msd_df %>%
  filter(MSD > 0) %>%  # skip tau=1 if noisy, avoid log(0)
  mutate(
    log_tau = log(tau),
    log_MSD = log(MSD)
  ) %>%
  group_by(temp, group, ID) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(log_MSD ~ log_tau, data = .x)),
    fit_summary = map(fit, tidy)
  ) %>%
  unnest(fit_summary) %>%
  filter(term == "log_tau") %>%
  select(temp, group, ID, estimate, std.error, p.value) %>%
  rename(
    alpha = estimate,
    alpha_se = std.error,
    alpha_p = p.value
  )
View(msd_slopes)
alpha_by_temp_group <- msd_slopes %>%
  #group_by(group,temp) %>%
  group_by(temp) %>%
  summarise(
    mean_alpha = mean(alpha, na.rm = TRUE),
    se_alpha = sd(alpha, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
View(alpha_by_temp_group)
write.csv(alpha_by_temp_group, "MSD_alpha_summary1.csv",row.names=TRUE)
write.csv(msd_slopes, "MSD_alpha_25V.csv",row.names=TRUE)

msd_by_temp_loglog <- msd_overall %>%
  filter(tau > 1, MSD_mean > 0) %>%
  group_by(temp, tau) %>%
  mutate(
    log_tau = log10(tau),
    log_MSD = log10(MSD_mean)
  )
msd_by_temp_loglog$temp <- as.factor(msd_by_temp_loglog$temp) 
p<-ggplot(msd_by_temp_loglog, aes(x = log_tau, y = log_MSD, color = as.factor(temp))) +
  geom_line(size = 1) +
  labs(
    title = "Average Log–Log MSD per Temperature",
    x = expression(log[10](tau)),
    y = expression(log[10](Mean~MSD)),
    color = "Temperature"
  ) +
  scale_color_manual(values = custom_colors) +  
  theme_minimal()
p


p<-ggplot(msd_slopes, aes(x = temp, y = alpha, group = temp)) +
  geom_boxplot(outlier.shape = NA, fill = "grey80", color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 1, color = "steelblue") +
  labs(
    x = "Temperature (°C)",
    y = expression(alpha),
    title = expression("Diffusion exponent"~alpha~"across temperatures")
  ) +
  theme_classic(base_size = 14)
p

#statistical test
msd_slopes$temp <- as.factor(msd_slopes$temp) 
alpha_tests <- msd_slopes %>%
  group_by(temp) %>%
  summarise(
    mean_alpha = mean(alpha, na.rm = TRUE),
    n = n(),
    t_test = list(t.test(alpha, mu = 1))
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    statistic = sapply(t_test, function(x) x$statistic)
  )

alpha_tests

alpha_tests <- msd_slopes %>%
  group_by(temp) %>%
  summarise(
    mean_alpha = mean(alpha, na.rm = TRUE),
    n = n(),
    t_above1 = list(t.test(alpha, mu = 1, alternative = "greater")),
    t_below1 = list(t.test(alpha, mu = 1, alternative = "less"))
  ) %>%
  mutate(
    p_above1 = sapply(t_above1, \(x) x$p.value),
    p_below1 = sapply(t_below1, \(x) x$p.value)
  )
alpha_tests

alpha_tests_flat <- alpha_tests %>%
  mutate(
    t_above1 = sapply(t_above1, function(x) paste(capture.output(print(x)), collapse = " ")),
    t_below1 = sapply(t_below1, function(x) paste(capture.output(print(x)), collapse = " "))
  )

write.csv(alpha_tests_flat, "alphatests.csv", row.names = FALSE)




View(msd_slopes)


# log-log plot ------------------------------------------------------------
# convert to log space for each point, remove invalid values
loglog_per_point <- msd_df %>%
  filter(MSD > 0, tau > 1) %>%  # avoid log(0) and tau=1 if noisy
  mutate(
    log_tau = log10(tau),
    log_MSD = log10(MSD)
  )

# Average log–log MSD per track, tau
loglog_per_track <- loglog_per_point %>%
  group_by(temp, ID, tau) %>%
  summarise(
    log_tau = mean(log_tau),  # constant across ID/tau anyway
    log_MSD = mean(log_MSD),
    .groups = "drop"
  )

#  average across tracks for each (temp, tau)
loglog_avg <- loglog_per_track %>%
  group_by(temp, tau, log_tau) %>%
  summarise(
    mean_log_MSD = mean(log_MSD, na.rm = TRUE),
    se_log_MSD = sd(log_MSD, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p <- ggplot(loglog_avg, aes(x = log_tau, y = mean_log_MSD, color = as.factor(temp))) +
  geom_line(size = 1) +
  geom_ribbon(aes(
    ymin = mean_log_MSD - se_log_MSD,
    ymax = mean_log_MSD + se_log_MSD,
    fill = as.factor(temp)
  ), alpha = 0.2, color = NA) +
  labs(
    title = "Mean Log–Log MSD Across Tracks",
    x = expression(log[10](tau)),
    y = expression(mean~log[10](MSD)),
    color = "Temperature",
    fill = "Temperature"
  ) +
  scale_color_viridis_e() +
  scale_fill_viridis_e() +
  theme_minimal()

p



# convex hull calculation-------------------------------------------------------------
library(purrr)
library(tibble)

convex_metrics <- df %>%
  group_by(group, temp, ID) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n > max_tau) %>%
  mutate(group_id = paste(group, ID, sep = "_")) %>%
  pmap_dfr(function(group, temp, ID, n, group_id) {
    # Extract track points
    track <- df %>%
      filter(group == !!group, ID == !!ID) %>%
      arrange(Frame)
    
    # if fewer than 3 points, area/perimeter = NA
    if(nrow(track) < 3) {
      return(tibble(
        group = group,
        temp = temp,
        ID = ID,
        convex_area = NA_real_,
        convex_perimeter = NA_real_
      ))
    }
    
    # convex hull
    hull_idx <- chull(track$Xpos, track$Ypos) #chull finds point indices that make it convex
    hull_points <- track[c(hull_idx, hull_idx[1]), ]  # close polygon-extract rows from track
    #have to do hull_idx[1] to close it and repeat that first point
    
    # shoelace formula computes polygon area from vertex coordinates
    #get all x except 1st one, all y except last one
    x <- hull_points$Xpos #make the vectors 
    y <- hull_points$Ypos
    area <- 0.5 * abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1]))
    
    # sum of distances between consecutive hull points 
    perimeter <- sum(sqrt(diff(x)^2 + diff(y)^2))
    
    tibble(
      group = group,
      temp = temp,
      ID = ID,
      convex_area = area,
      convex_perimeter = perimeter
    )
  })
View(convex_metrics)
write.csv(convex_metrics, "raw_convex_met.csv",row.names=TRUE)
convex_summary <- convex_metrics %>%
  group_by(temp) %>%
  summarise(
    mean_area = mean(convex_area, na.rm = TRUE),
    sd_area = sd(convex_area, na.rm = TRUE),
    se_area = sd_area / sqrt(sum(!is.na(convex_area))),
    
    median_area = median(convex_area, na.rm = TRUE),
    iqr_area = IQR(convex_area, na.rm = TRUE),
    
    mean_perimeter = mean(convex_perimeter, na.rm = TRUE),
    sd_perimeter = sd(convex_perimeter, na.rm = TRUE),
    se_perimeter = sd_perimeter / sqrt(sum(!is.na(convex_perimeter))),
    
    median_perimeter = median(convex_perimeter, na.rm = TRUE),
    iqr_perimeter = IQR(convex_perimeter, na.rm = TRUE),
    
    n_tracks = n(),
    .groups = "drop"
  )
View(convex_summary)
write.csv(convex_summary, "convex_metrics.csv",row.names=TRUE)


# convex hull binned by AR ------------------------------------------------
# 1. Add shape bins to your df (motility_mode2)
# Step 1: Bin raw AR values

library(purrr)
library(tibble)


# 1. Add AR bins to per-frame data
df_binned <- df %>%
  mutate(
    AR_bin = case_when(
      AR >= 1   & AR < 1.1   ~ "1-1.1",
      AR >= 1.1 & AR < 1.25  ~ "1.1-1.25",
      AR >= 1.25 & AR < 1.5  ~ "1.25-1.5",
      AR >= 1.5 & AR < 2.0   ~ "1.5-2",
      AR >= 2.0 & AR < 2.5   ~ "2-2.5",
      AR >= 2.5 & AR < 3     ~ "2.5-3",
      AR >= 3                ~ "3.0",
      TRUE ~ "Other"
    )
  )

epsilon <- 1e-8

# mean_AR per track
track_summary <- df %>%
  group_by(temp, group, ID) %>%
  summarise(
    mean_AR = mean(AR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    AR_shape = case_when(
      mean_AR >= 1 & mean_AR < 1.25 & temp %in% c(25, 70)          ~ "cyst",
      mean_AR >= 1.25 & mean_AR < 2.5 & temp %in% c(55,57, 60, 62, 64)  ~ "lamelli",
      mean_AR >= 2.5 & temp %in% c(55,57, 60, 62, 64)                   ~ "worm",
      TRUE ~ NA_character_
    )
  )

# convex hull metrics — per whole track
convex_metrics <- track_summary %>%
  pmap_dfr(function(temp, group, ID, mean_AR, AR_shape) {
    track <- df %>%
      filter(temp == !!temp, group == !!group, ID == !!ID, stepwise_displacement < 0.55) %>%
      arrange(Frame)
    
    if(nrow(track) < 3) {
      return(tibble(
        group = group,
        temp = temp,
        ID = ID,
        AR_shape = AR_shape,
        convex_area = NA_real_,
        convex_perimeter = NA_real_
      ))
    }
    
    # convex hull
    hull_idx <- chull(track$Xpos, track$Ypos)
    hull_points <- track[c(hull_idx, hull_idx[1]), ]
    x <- hull_points$Xpos
    y <- hull_points$Ypos
    area <- 0.5 * abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1]))
    perimeter <- sum(sqrt(diff(x)^2 + diff(y)^2))
    
    tibble(
      group = group,
      temp = temp,
      ID = ID,
      AR_shape = AR_shape,
      convex_area = area,
      convex_perimeter = perimeter
    )
  })

# summarize by temp + shape
convex_summary <- convex_metrics %>%
  group_by(temp, AR_shape) %>%
  summarise(
    mean_area = mean(convex_area, na.rm = TRUE),
    sd_area = sd(convex_area, na.rm = TRUE),
    se_area = sd_area / sqrt(sum(!is.na(convex_area))),
    median_area = median(convex_area, na.rm = TRUE),
    iqr_area = IQR(convex_area, na.rm = TRUE),
    mean_perimeter = mean(convex_perimeter, na.rm = TRUE),
    sd_perimeter = sd(convex_perimeter, na.rm = TRUE),
    se_perimeter = sd_perimeter / sqrt(sum(!is.na(convex_perimeter))),
    median_perimeter = median(convex_perimeter, na.rm = TRUE),
    iqr_perimeter = IQR(convex_perimeter, na.rm = TRUE),
    n_tracks = n(),
    .groups = "drop"
  )

View(convex_summary)
write.csv(convex_summary,"convexhull_byAR.csv")

