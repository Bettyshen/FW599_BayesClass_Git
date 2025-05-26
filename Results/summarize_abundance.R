# Goal: summarize the abundance results
# Author: Fang-Yu (Betty) Shen
# Date: 2025-05-26

# Load required packages
library(dplyr)
library(tidyr)
library(knitr)
library(coda)

# Read MCMC output if not already in environment
mcmc.output <- readRDS("/Users/shenf/Documents/Bayes_Class/Final_project/Results/mcmc_output.rds")


#Calculate proportion of posterior draw greater or lesser than zero
get_parameter_summary <- function(samples_matrix) {
  summary_list <- apply(samples_matrix, 2, function(x) {
    median_val <- median(x)
    cri <- unname(quantile(x, c(0.025, 0.975)))  # remove names
    prob_pos <- mean(x > 0)
    prob_neg <- mean(x < 0)
    prob <- max(prob_pos, prob_neg)
    direction <- ifelse(prob_pos > prob_neg, "positive", "negative")
    c(median = median_val,
      lower95 = cri[1],
      upper95 = cri[2],
      prob_diff_from_zero = prob,
      direction = direction)
  })

  df <- as.data.frame(t(summary_list))
  df$median <- as.numeric(df$median)
  df$lower95 <- as.numeric(df$lower95)
  df$upper95 <- as.numeric(df$upper95)
  df$prob_diff_from_zero <- as.numeric(df$prob_diff_from_zero)
  df$direction <- as.character(df$direction)
  return(df)
}


# Convert MCMC output to matrix format
if(inherits(mcmc.output, "mcmc.list")) {
  N_samples <- do.call(rbind, mcmc.output)
} else if(inherits(mcmc.output$samples, "mcmc.list")) {
  N_samples <- do.call(rbind, mcmc.output$samples)
} else {
  N_samples <- as.matrix(mcmc.output$samples)
}

# Get N parameters and summarize
N_cols <- grep("^N\\[", colnames(N_samples), value = TRUE)
N_matrix <- N_samples[, N_cols]

# Get summary for all N parameters
N_summary <- get_parameter_summary(N_matrix)

# Calculate overall statistics
overall_summary <- data.frame(
  Statistic = c(
    "Total Population (sum of medians)",
    "Mean abundance per site",
    "Median abundance per site",
    "Sites with high probability (>0.95) of at least one individual",
    "Minimum site abundance",
    "Maximum site abundance"
  ),
  Value = c(
    sum(N_summary$median),
    mean(N_summary$median),
    median(N_summary$median),
    sum(N_summary$prob_diff_from_zero > 0.95),
    min(N_summary$median),
    max(N_summary$median)
  )
)
# Write overall summary to csv
write.csv(overall_summary, "/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_summary.csv", row.names = FALSE)

# Calculate abundance categories
abundance_categories <- cut(N_summary$median,
                          breaks = c(0, 1, 5, 10, 20, Inf),
                          labels = c("0-1", "1-5", "5-10", "10-20", ">20"))
category_summary <- data.frame(
  Category = levels(abundance_categories),
  Count = as.numeric(table(abundance_categories)),
  Proportion = as.numeric(prop.table(table(abundance_categories)))
)

# Print results
cat("\nOverall Summary Statistics:\n")
print(knitr::kable(overall_summary, digits = 2))

cat("\nDistribution of Site Abundances:\n")
print(knitr::kable(category_summary, digits = 3))

# Calculate quantiles of site abundances
abundance_quantiles <- data.frame(
  Quantile = c("2.5%", "25%", "50%", "75%", "97.5%"),
  Value = quantile(N_summary$median, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
)

cat("\nQuantiles of Site Abundances:\n")
print(knitr::kable(abundance_quantiles, digits = 2))

# Save results
results <- list(
  N_summary = N_summary,
  overall_summary = overall_summary,
  category_summary = category_summary,
  abundance_quantiles = abundance_quantiles
)


saveRDS(results, 
        file = "/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_summary.rds")

# Optional: Create a histogram of median abundances
if(require(ggplot2)) {
  p <- ggplot(data.frame(median = N_summary$median), aes(x = median)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(x = "Median Abundance", y = "Number of Sites",
         title = "Distribution of Estimated Site Abundances") +
    theme_bw()
  
  ggsave("/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_histogram.png", 
         p, width = 8, height = 6, dpi = 300)
} 

#===============================================
# Create comparison plots
if(require(ggplot2)) {
  # Prepare data for plotting
  # Get eBird raw counts
  raw_eBird <- sapply(1:nsite, function(i) {
    max(eBird_cnt[which(sites_eBird == i)], na.rm = TRUE)
  })
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Abundance = c(N_summary$median, raw_eBird),
    Type = factor(rep(c("Model Estimated", "eBird Raw"), 
                     c(length(N_summary$median), length(raw_eBird))),
                 levels = c("Model Estimated", "eBird Raw"))
  )
  
  # Create histogram
  p1 <- ggplot(plot_data, aes(x = Abundance, fill = Type)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
    scale_fill_manual(values = c("Model Estimated" = "skyblue", 
                                "eBird Raw" = "darkgreen")) +
    scale_x_continuous(breaks = 0:11, limits = c(0, 11)) +
    labs(x = "Abundance", y = "Number of Sites",
         title = "Distribution of eBird Raw Abundance vs. Estimated Abundance") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(size = 12),
          axis.text.x = element_text(angle = 0))
  
  # Create density plot
  p2 <- ggplot(plot_data, aes(x = Abundance, fill = Type)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("Model Estimated" = "skyblue", 
                                "eBird Raw" = "darkgreen")) +
    scale_x_continuous(breaks = 0:11, limits = c(0, 11)) +
    labs(x = "Abundance", y = "Density",
         title = "Density Distribution of eBird Raw Abundance vs. Estimated Abundance") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(size = 12),
          axis.text.x = element_text(angle = 0))

  # Save plots
  ggsave("/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_comparison_histogram.png", 
         p1, width = 10, height = 7, dpi = 300)
  ggsave("/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_comparison_density.png", 
         p2, width = 10, height = 7, dpi = 300)
  
  # Calculate summary statistics for comparison
  comparison_summary <- data.frame(
    Metric = c("Mean", "Median", "Maximum", "Proportion Zeros",
               "25th Percentile", "75th Percentile"),
    Model_Estimated = c(mean(N_summary$median), 
                       median(N_summary$median), 
                       max(N_summary$median),
                       mean(N_summary$median == 0),
                       quantile(N_summary$median, 0.25),
                       quantile(N_summary$median, 0.75)),
    eBird_Raw = c(mean(raw_eBird), 
                  median(raw_eBird), 
                  max(raw_eBird),
                  mean(raw_eBird == 0),
                  quantile(raw_eBird, 0.25),
                  quantile(raw_eBird, 0.75))
  )
  
  # Write comparison summary to csv
  write.csv(comparison_summary, 
            "/Users/shenf/Documents/Bayes_Class/Final_project/Results/abundance_comparison_summary.csv", 
            row.names = FALSE)
  
  # Print comparison summary
  cat("\nComparison of eBird Raw Counts vs. Model-Estimated Abundances:\n")
  print(knitr::kable(comparison_summary, digits = 3))
} 
