#===============================================
# Integrated population model (Oregon 2020 + eBird) for California Scrub-Jay
# Reference: Zhao, Q., Latif, Q. S., Nuse, B. L., Pavlacky Jr, D. C., Kilner, C. L., Ryder, T. B., & Latimer, C. E. (2024).
# Integrating counts from rigorous surveys and participatory science to better understand spatiotemporal variation in population processes. 
# Methods in Ecology and Evolution, 15(8), 1380-1393.

#===============================================
rm(list=ls())
library(nimble)
library(rstan)
library(data.table)
library(dplyr)
library(ggplot2)
library(coda)
library(MCMCvis)
library(tidyverse)
library(mcmcplots)
library(ragg)


setwd("/Users/shenf/Documents/Bayes_Class/Final_project/")
source('/Users/shenf/Documents/Bayes_Class/Week2_Nimble/attach.nimble_v2.R')
#=============================================
# Read in data
#=============================================

# Scrub-Jay data
OR2020 <- read.csv("/Users/shenf/Documents/Bayes_Class/Final_project/Data/CaliScrubJay_OR2020.csv")
eBird <- read.csv("/Users/shenf/Documents/Bayes_Class/Final_project/Data/CaliScrubJay_eBird.csv")
OR2020_cnt <- read.csv("/Users/shenf/Documents/Bayes_Class/Final_project/Data/OR2020_cnt.csv")

#=============================================
# Organize data
#=============================================

# ===== OR2020 data ======= #
nsite <- dim(OR2020)[1] # Number of OR2020 sites (1862)
ndist <- 20 # Number of distance bins
ntime <- 6 # Number of time bins
cutoff <- 200 # Maximum distance for distance sampling & future density estimates
breaks <- seq(0, cutoff, length.out=ndist + 1) # total distance bins for calculating perceptibility
nobs_OR2020 <- dim(OR2020_cnt)[1] # number of observations

# Create new sequential IDs for OR2020
OR2020$new_id <- 1:nsite
OR2020_cnt$new_id <- 1:nobs_OR2020

# Create site mapping
info <- data.frame(site = 1:nsite)

# Create site indices that align with observations
sites_OR2020 <- OR2020_cnt$new_id

# Fix OR2020_cnt structure (observations x time interval + distance bins)
# Remove the ID column and convert to numeric matrix
OR2020_cnt <- as.matrix(OR2020_cnt[,-1])  # Remove the first column from old ID
OR2020_cnt <- matrix(as.numeric(OR2020_cnt), nrow=nobs_OR2020, ncol=ndist*ntime)  # Convert to numeric

# Create distance-time interval names
# First create all combinations
dist_bins <- seq(0, cutoff-10, by=10)
time_bins <- 1:ntime

# Create the names in the exact order (time interval x distance bins)
dist_time_names <- character(ndist * ntime)
idx <- 1
for(t in 1:ntime) {
  for(d in 1:ndist) {
    dist_time_names[idx] <- paste0("X", dist_bins[d], ".", dist_bins[d]+10, "..", t)
    idx <- idx + 1
  }
}

# Create observation names (using site indices)
obs_names <- as.character(1:nobs_OR2020)

# Add dimnames to the matrix
dimnames(OR2020_cnt) <- list(obs_names, dist_time_names)

# Calculate the total number of observations for each site
OR2020_cnt_sum <- rowSums(OR2020_cnt)

# ========= eBird data ==============#
# create new sequential IDs
eBird <- eBird[sample(1:dim(eBird)[1], nsite),]
eBird$new_id <- 1:dim(eBird)[1]
nobs_eBird <- dim(eBird)[1]
sites_eBird <- eBird$new_id
locat_eBird <- sites_eBird  # Use the same mapping for locations
nlocat_eBird <- length(unique(locat_eBird))
duration_eBird <- log(eBird$duration_minutes)
eBird_cnt <- eBird$observation_count

# ======== Environmental data =========#
precip <- OR2020$ppt_median
EVI <- OR2020$EVI_median
tmean <- OR2020$tmean_median
min.OR2020 <- OR2020$MinutesAfterSunrise
day.OR2020 <- OR2020$Day

# ======== eBird covariates ==========#
# For perceptibility & availability calculation
min.eBird <- eBird$MinutesAfterSunrise
day.eBird <- eBird$Day

# Print data summaries to check for issues
cat("Number of OR2020 sites:", nsite, "\n")
cat("Number of OR2020 observations:", nobs_OR2020, "\n")
cat("Number of eBird observations:", nobs_eBird, "\n")
cat("Number of unique eBird locations:", nlocat_eBird, "\n")
cat("Range of eBird counts:", range(eBird_cnt, na.rm = TRUE), "\n")
cat("Range of duration:", range(duration_eBird, na.rm = TRUE), "\n")
cat("Range of minutes after sunrise:", range(min.eBird, na.rm = TRUE), "\n")

#=============================================
# Define model in Nimble for Scrub-Jay data 
#=============================================
code <- nimbleCode({

  # Priors
  zeta ~ dunif(0, 1) # Occupancy
  for (k in 1:4) {
    beta_lambda0[k] ~ dnorm(0, sd=10) # Bird-location covariates
  } # k
  sigma ~ dgamma(.01, .01) # OR2020 half-normal distance function shape parameter
  theta ~ dunif(0, 1) # Availability within one time interval

  for (k in 1:4) {
    omega_beta[k] ~ dnorm(0, sd=10) # Effects of bird-location covariates on eBird availability & observer perceptibility
  } # k
  omega_zeta ~ dunif(0, 1) # eBird Occupancy
  omega_delta ~ dgamma(.01, .01) # eBird Noise
  omega_locat_sd ~ dgamma(0.01, .01) # location SD

  # ===== Process model ===== #
  for (i in 1:nsite) {
    z[i] ~ dbinom(zeta, 1)
    lambda0[i] <- exp(
      beta_lambda0[1] + 
      beta_lambda0[2] * precip[i] + 
      beta_lambda0[3] * EVI[i] + 
      beta_lambda0[4] * tmean[i])
    N[i] ~ dpois(lambda0[i] * z[i] + 5e-3 * (1 - z[i]))
  } # i

  # ===== Observation model (likelihood) - Oregon 2020 ===== #
  for (d in 1:ndist) {
    intval[d] <- sigma^2 * (exp(-1*(breaks[d]^2)/(2*sigma^2)) - exp(-1*(breaks[d+1]^2)/(2*sigma^2)))
    psi[d] <- 2 * intval[d] / (cutoff ^ 2) # Perceptibility
  } # d
  psi_sum <- sum(psi[1:ndist]) # Sum of perceptibility across all distance bins
  for(d in 1:ndist) {
    psi_prop[d] <- psi[d] / psi_sum # Proportion of perceptibility for each distance bin
  } # d

  for (r in 1:ntime) {
    phi[r] <- (1 - theta) ^ (r - 1) * theta # Availability
  } # r
  phi_sum <- sum(phi[1:ntime]) # Sum of availability across all time intervals
  for (r in 1:ntime) {
    phi_prop[r] <- phi[r] / phi_sum # Proportion of availability for each time interval
  } # r

  # Calculate pi at specific distance bin and time interval (detection probability = perceptibility * availability)
  for (d in 1:ndist) {
    for (r in 1:ntime) {
      pi_raw[(r-1)*ndist+d] <- psi_prop[d] * phi_prop[r] # Detection probability
    } # r
  } # d

  # Add small constant and normalize
  pi_sum <- sum(pi_raw[1:(ndist*ntime)])
  for(k in 1:(ndist*ntime)) {
    pi[k] <- (pi_raw[k] + 1e-10) / (pi_sum + (ndist*ntime)*1e-10) # Normalized detection probability
  }
  # Expected count at specific distance bin and time interval
  for(k in 1:nobs_OR2020) {
    OR2020_cnt_sum[k] ~ dbinom(psi_sum * phi_sum, N[sites_OR2020[k]])
    OR2020_cnt[k,1:(ndist*ntime)] ~ dmultinom(pi[1:(ndist*ntime)], OR2020_cnt_sum[k])
  } # k

  # ===== Observation model (likelihood) - eBird ===== #
  for (k in 1:nlocat_eBird) {
    omega_locat_epsilon[k] ~ dnorm(0, sd = omega_locat_sd) # Random effect for location
  } #k
  for (k in 1:nobs_eBird) {
    omega_z[k] ~ dbinom(omega_zeta, 1) # Occupancy
    omega[k] <- exp(
      omega_beta[1] + 
      omega_beta[2] * EVI[sites_eBird[k]] + 
      omega_beta[3] * min.eBird[k] +
      omega_beta[4] * duration_eBird[k] +
      omega_locat_epsilon[locat_eBird[k]] # Random effect for location
    )

    # Expected count from eBird
    eBird_cnt[k] ~ dpois(N[sites_eBird[k]] * omega[k] * omega_z[k] + omega_delta * (1 - omega_z[k])) 
  } # i

}) # nimbleCode

#============  Run Nimble model  =============#
# Scrub-Jay data
constants <- list(
  nsite=nsite, 
  ndist=ndist, ntime=ntime, breaks=breaks, cutoff=cutoff, 
  nobs_OR2020=nobs_OR2020, sites_OR2020=sites_OR2020, 
  nobs_eBird=nobs_eBird, sites_eBird=sites_eBird, 
  nlocat_eBird=nlocat_eBird, locat_eBird=locat_eBird, 
  duration_eBird=c(scale(duration_eBird)),
  precip=c(scale(precip)),
  EVI=c(scale(EVI)),
  tmean=c(scale(tmean)),
  min.eBird=c(scale(min.eBird))
)

data <- list(
  OR2020_cnt_sum=OR2020_cnt_sum, OR2020_cnt=OR2020_cnt, 
  eBird_cnt=eBird_cnt
)

# MCMC Settings
ni <- 100000
nt <- 10
nb <- 50000
nc <- 3

parameters <- c("zeta", "beta_lambda0", "sigma", "theta", 
"omega_beta", "omega_zeta", "omega_delta", "omega_locat_sd", "N")

# Initial values
Ni <- rep(5, nsite)  # Changed from matrix to vector
for (i in 1:nsite) {
  if (length(c(OR2020_cnt_sum[which(sites_OR2020 == i)], eBird_cnt[which(sites_eBird == i)])) > 0) {
    Ni[i] <- max(c(OR2020_cnt_sum[which(sites_OR2020 == i)], eBird_cnt[which(sites_eBird == i)]))
  }
} # i

zi <- ifelse(Ni > 0, 1, 0)  # Changed from rowSums to direct comparison
Ni <- Ni * 2 + 10

# Create initial values for the random effect
omega_locat_epsilon_init <- rnorm(nlocat_eBird, 0, 0.1)

inits <- list(
  zeta = 0.2, 
  beta_lambda0 = c(2, 0.5, -1, -1),  # Adjusted for scaled covariates
  sigma = 60, 
  theta = 0.4, 
  omega_beta = c(-6, 0.1, 0.5, 0.2),  # Adjusted for scaled covariates
  omega_zeta = 0.3,
  omega_z = ifelse(eBird_cnt>0, 1, 0),
  omega_delta = 0.01,
  omega_locat_sd = 0.1,
  omega_locat_epsilon = omega_locat_epsilon_init,
  N = Ni,  
  z = zi
)

mcmc.output <- nimbleMCMC(code = code,
                         data = data,
                         constants = constants,
                         inits = inits,  
                         monitors = parameters,
                         niter = ni,
                         nburnin = nb,
                         nchains = nc,
                         thin = nt,
                         summary = TRUE,
                         samplesAsCodaMCMC = TRUE)

# Save MCMC output
saveRDS(mcmc.output, file = "/Users/shenf/Documents/Bayes_Class/Final_project/Results/mcmc_output.rds")

# MCMC plot

attach.nimble(mcmc.output$samples)
mcmcplot(mcmc.output$samples)

# Examine MCMC output
# 95% Credible interval
  # N
quantile(N, c(0.025, 0.975))  # Overall population size across all sites
apply(N, 2, quantile, c(0.025, 0.5, 0.975))  # Site-specific estimates


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

# Apply it to all MCMC samples
summary_df <- get_parameter_summary(as.matrix(mcmc.output$samples))
print(summary_df)

# 95% CRI for parameters
write.csv(summary_df, "/Users/shenf/Documents/Bayes_Class/Final_project/Results/parameter_summary.csv", row.names = FASLE)
