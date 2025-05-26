# Goal: plot the posterior distributions of the parameters
# Author: Fang-Yu (Betty) Shen
# Date: 2025-05-26

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(coda)

# Read MCMC output
mcmc.output <- readRDS("/Users/shenf/Documents/Bayes_Class/Final_project/Results/mcmc_output.rds")

# Get all parameter names except N
param_names <- c("zeta", "beta_lambda0", "sigma", "theta", 
                "omega_beta", "omega_zeta", "omega_delta", "omega_locat_sd")

# Convert MCMC output to matrix format
if(inherits(mcmc.output, "mcmc.list")) {
  samples <- do.call(rbind, mcmc.output)
} else if(inherits(mcmc.output$samples, "mcmc.list")) {
  samples <- do.call(rbind, mcmc.output$samples)
} else {
  samples <- as.matrix(mcmc.output$samples)
}

# Get columns that match our parameters (excluding N)
param_cols <- grep(paste(paste0("^", param_names, "\\[?"), collapse="|"), 
                  colnames(samples), value=TRUE)

# Create a data frame for plotting
plot_data <- as.data.frame(samples[, param_cols])

# Reshape data for plotting
plot_data_long <- plot_data %>%
  tidyr::pivot_longer(cols = everything(),
                     names_to = "Parameter",
                     values_to = "Value")

# Create more readable parameter names
plot_data_long <- plot_data_long %>%
  mutate(Parameter = case_when(
    grepl("^beta_lambda0\\[1\\]", Parameter) ~ "β₀ (Intercept)",
    grepl("^beta_lambda0\\[2\\]", Parameter) ~ "β₁ (Precipitation)",
    grepl("^beta_lambda0\\[3\\]", Parameter) ~ "β₂ (Enhanced Vegetation Index)",
    grepl("^beta_lambda0\\[4\\]", Parameter) ~ "β₃ (Temperature)",
    grepl("^omega_beta\\[1\\]", Parameter) ~ "ω₀ (eBird Intercept)",
    grepl("^omega_beta\\[2\\]", Parameter) ~ "ω₁ (eBird Enhanced Vegetation Index)",
    grepl("^omega_beta\\[3\\]", Parameter) ~ "ω₂ (eBird Minutes since sunrise)",
    grepl("^omega_beta\\[4\\]", Parameter) ~ "ω₃ (eBird Duration)",
    grepl("^zeta", Parameter) ~ "ζ (Occupancy)",
    grepl("^sigma", Parameter) ~ "σ (Detection function shape parameter)",
    grepl("^theta", Parameter) ~ "θ (Availability within one time interval)",
    grepl("^omega_zeta", Parameter) ~ "ωζ (eBird Occupancy)",
    grepl("^omega_delta", Parameter) ~ "ωδ (eBird Noise)",
    grepl("^omega_locat_sd", Parameter) ~ "ωσ (Location SD)",
    TRUE ~ Parameter
  ))

# Create faceted density plot
p <- ggplot(plot_data_long, aes(x = Value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  facet_wrap(~Parameter, scales = "free", ncol = 3) +
  labs(x = "Parameter Value",
       y = "Density",
       title = "Posterior Distributions of Parameters") +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(size = 12, hjust = 0.5))

# Save the plot
ggsave("/Users/shenf/Documents/Bayes_Class/Final_project/Results/parameter_posteriors.png", 
       p, width = 12, height = 8, dpi = 300)

# Select parameters to plot (for presentation) =================================
# Create more readable parameter names

# Reshape data for plotting
plot_data_selected <- plot_data %>%
  tidyr::pivot_longer(cols = everything(),
                     names_to = "Parameter",
                     values_to = "Value") %>%
  filter(grepl("^beta_lambda0\\[[234]\\]|^omega_beta\\[[234]\\]", Parameter)) %>%  # Only keep selected parameters
  mutate(Parameter = case_when(
    grepl("^beta_lambda0\\[2\\]", Parameter) ~ "β₁ (Precipitation)",
    grepl("^beta_lambda0\\[3\\]", Parameter) ~ "β₂ (Enhanced Vegetation Index)",
    grepl("^beta_lambda0\\[4\\]", Parameter) ~ "β₃ (Temperature)",
    grepl("^omega_beta\\[2\\]", Parameter) ~ "ω₁ (eBird Enhanced Vegetation Index)",
    grepl("^omega_beta\\[3\\]", Parameter) ~ "ω₂ (eBird Minutes since sunrise)",
    grepl("^omega_beta\\[4\\]", Parameter) ~ "ω₃ (eBird Duration)",
    TRUE ~ Parameter
  ))

# Create faceted density plot
p.selected <- ggplot(plot_data_selected, aes(x = Value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  facet_wrap(~Parameter, scales = "free", ncol = 2) +
  labs(x = "Parameter Value",
       y = "Density",
       title = "Posterior Distributions of Selected Parameters") +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(size = 14, hjust = 0.5))

# Save the plot
ggsave("/Users/shenf/Documents/Bayes_Class/Final_project/Results/parameter_posteriors_selected.png", 
       p.selected, width = 10, height = 8, dpi = 300)
