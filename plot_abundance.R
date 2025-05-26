# Goal: plot the estimated abundance vs. raw counts
# Author: Fang-Yu (Betty) Shen
# Date: 2025-05-26

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)

# Read the MCMC output if not already in environment
mcmc.output <- readRDS("/Users/shenf/Documents/Bayes_Class/Final_project/Results/mcmc_output.rds")

# First, let's examine the structure of the MCMC output
print("Structure of MCMC output:")
str(mcmc.output)

print("Column names in MCMC samples:")
print(colnames(mcmc.output$samples))

# Get the N parameters
N_cols <- grep("^N\\[", colnames(mcmc.output$samples), value = TRUE)
print("N parameter columns found:")
print(N_cols)

# Calculate median and 95% CI for N at each site
N_matrix <- as.matrix(mcmc.output$samples[, N_cols])
N_estimates <- apply(N_matrix, 2, 
                    function(x) c(median = median(x), 
                                lower = quantile(x, 0.025),
                                upper = quantile(x, 0.975)))
N_estimates <- t(N_estimates)

# Create a data frame with the estimates
results_df <- data.frame(
  site = 1:nrow(N_estimates),
  estimated_N = N_estimates[, "median"],
  lower_CI = N_estimates[, "lower"],
  upper_CI = N_estimates[, "upper"]
)

# Add raw counts
results_df$OR2020_raw <- sapply(1:nsite, function(i) {
  max(OR2020_cnt_sum[which(sites_OR2020 == i)], na.rm = TRUE)
})
results_df$eBird_raw <- sapply(1:nsite, function(i) {
  max(eBird_cnt[which(sites_eBird == i)], na.rm = TRUE)
})

# Create long format for plotting
plot_df <- results_df %>%
  pivot_longer(cols = c(OR2020_raw, eBird_raw),
               names_to = "data_source",
               values_to = "raw_count")

# Create the plot
p <- ggplot(plot_df, aes(x = raw_count, y = estimated_N)) +
  geom_point(aes(color = data_source), alpha = 0.5) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("OR2020_raw" = "blue", "eBird_raw" = "red"),
                    labels = c("OR2020", "eBird")) +
  labs(x = "Raw Count",
       y = "Estimated Abundance (N)",
       color = "Data Source",
       title = "Estimated Abundance vs. Raw Counts",
       subtitle = "Dashed line represents 1:1 relationship") +
  theme_bw() +
  theme(legend.position = "bottom")

# Save the plot
ggsave("/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_comparison.png", p, width = 10, height = 8, dpi = 300)

# Print summary statistics
cat("\nSummary of Estimates vs. Raw Counts:\n")
cat("Median estimated N:", median(results_df$estimated_N), "\n")
cat("Median OR2020 count:", median(results_df$OR2020_raw), "\n")
cat("Median eBird count:", median(results_df$eBird_raw), "\n") 