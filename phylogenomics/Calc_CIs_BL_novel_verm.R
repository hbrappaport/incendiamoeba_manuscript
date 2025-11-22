library(ggplot2)

## branch length comparison

bl = read.csv('phylogenomics/branch_length_comparison.csv')

# group data
group = bl$Avg_CA[1:6]

# obs value to test
obs = bl$Avg_CA[7]

# t-test
t_test_result = t.test(group, mu = obs)
t_test_result

# observed value lie within the 95% confidence interval?
ci = t_test_result$conf.int
cat("Obs value within 95% CI:", obs >= ci[1] && obs <= ci[2], "\n")

# visualize
df = data.frame(value = group)

# ggplot
ggplot(df, aes(x = value)) +
  geom_histogram(binwidth = 0.012, fill = "gray", color = "gray", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(value)), color = "black", size = 0.5) +
  geom_vline(xintercept = obs, color = "red", size = 0.5) +
  geom_vline(xintercept = ci[1], color = "#666666", linetype = "dashed", size = 0.25) +
  geom_vline(xintercept = ci[2], color = "#666666", linetype = "dashed", size = 0.25) +
  labs(x = "Value", y = "Count") +
  theme_classic()
