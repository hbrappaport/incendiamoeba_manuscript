library(dplyr)
library(viridis) #color aesthetics
library(ggbeeswarm) 
library(tidyr)
library(readr)

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

df<-data
#df$temp <- "66"  # name group in df

summary_stats <- df %>%
  group_by(temp, group,ID) %>%  # group by temp group and ID to get particle track
  filter(stepwise_speed > 0, stepwise_speed < 10) %>%  # quality control for detached cells
  summarise(
    mean_speed = (mean(stepwise_speed, na.rm = TRUE))*60, #per minute
    mean_area = (mean(area, na.rm = TRUE)),
    median_area = (median(area, na.rm = TRUE)),
    median_speed=median(stepwise_speed,na.rm=TRUE)*60,
    sd_speed = sd(stepwise_speed, na.rm = TRUE),
    mean_AR = mean(AR, na.rm = TRUE),
    median_AR=median(AR,na.rm=TRUE),
    sd_AR = sd(AR, na.rm = TRUE),
    mean_major=mean(Major,na.rm=TRUE),
    mean_minor=mean(Minor,na.rm=TRUE),
    mean_PARatio = mean(PARatio, na.rm = TRUE),
    sd_PARatio = sd(PARatio, na.rm = TRUE),
    distance=(mean(stepwise_speed,na.rm=TRUE))*60, #per minute
    n = n()  # Count of observations
  ) %>%
  ungroup() %>%
  mutate(
    mean_speed = mean_speed,
    ci_speed_upper = mean_speed + (1.645 * sd_speed / sqrt(n)),
    ci_speed_lower = mean_speed - (1.645 * sd_speed / sqrt(n))
  ) %>%
  filter(mean_speed >= ci_speed_lower & mean_speed <= ci_speed_upper) 

IQRstats <- summary_stats %>%
  group_by(temp) %>% 
  summarise(
    median_speed2=median(median_speed,na.rm=TRUE),
    speed_IQR = IQR(median_speed, na.rm = TRUE),
    speed_Q1 = quantile(median_speed, 0.25, na.rm = TRUE),
    speed_Q3 = quantile(median_speed, 0.75, na.rm = TRUE),
    speed_p5 = quantile(median_speed, 0.05, na.rm = TRUE),
    speed_p95 = quantile(median_speed, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()
View(IQRstats)

write.csv(IQRstats,"IQR_Velocity_stats_withVerm.csv")

View(summary_stats)
mean_stats <- summary_stats %>%
  group_by(temp) %>%  # group by temp, this should already be reported per track 
  summarise(
    mean_speed = (mean(mean_speed, na.rm = TRUE)),
    mean_area = (mean(mean_area, na.rm = TRUE)),
    median_area = (median(median_area, na.rm = TRUE)),
    median_speed=median(median_speed,na.rm=TRUE),
    mean_major=mean(mean_major,na.rm=TRUE),
    mean_minor=mean(mean_minor,na.rm=TRUE),
    mean_AR = mean(mean_AR, na.rm = TRUE),
    median_AR=median(median_AR,na.rm=TRUE),
    n = n()  # count of observations
  ) %>%
  ungroup() 
View(mean_stats)

#summary_stats$temp <- "57"  # name group in df
write.csv(mean_stats,"mean_stats.csv")
#clear console -> cat("\014")

# Visualization -----------------------------------------------------------
summary_stats$temp <- as.factor(summary_stats$temp)  #  convert to factor
p<- ggplot(summary_stats, aes(x=temp, y=mean_speed,group=temp))+ coord_cartesian(ylim=c(0,95)) + geom_boxplot(outlier.shape = NA)+ 
  geom_beeswarm(dodge.width=1,aes(color = temp), size = 0.5, alpha = 0.5,cex=0.25) + 
  #scale_color_viridis_d(option = "D", direction = 1)#color by group eventully
  scale_color_manual(values = custom_colors)
p 
p<-p + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank()) + 
  ylab("Median Velocity µm/min") 
p
