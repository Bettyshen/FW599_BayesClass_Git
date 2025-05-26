# Goal: transform Oregon 2020 data to match Qing (2024) format
#==============
# Read in data
#==============
# Scrub-Jay data
OR2020 <- read.csv("/Users/shenf/Documents/Bayes_Class/Final_project/Data/CaliScrubJay_OR2020.csv")

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# Define breakpoints
dist_breaks <- seq(0, 200, by = 10)
time_bins <- 1:6

# Create column names for output
bin_labels <- expand.grid(
  dist = paste0(dist_breaks[-length(dist_breaks)], ".", dist_breaks[-1]),
  time = time_bins
) %>%
  mutate(colname = paste0("X", dist, "..", time)) %>%
  pull(colname)

# Initialize output matrix
out_mat <- matrix(0, nrow = nrow(OR2020), ncol = length(bin_labels))
colnames(out_mat) <- bin_labels

# Function to find the first interval with detection
get_first_detection_interval <- function(x) {
  int_idx <- which(x == 1)
  if (length(int_idx) > 0) return(int_idx[1]) else return(NA_integer_)
}

# Loop over each row of the data
for (i in seq_len(nrow(OR2020))) {
  det_time <- get_first_detection_interval(as.numeric(OR2020[i, grep("^Interval_", names(OR2020))]))
  dist <- OR2020$Distance[i]
  
  # Skip if no time or no distance info
  if (is.na(det_time) || is.na(dist)) next
  
  # Assign to the correct bin
  dist_bin <- findInterval(dist, dist_breaks, rightmost.closed = TRUE)
  
  # Don't assign if outside defined bins (e.g., >120 m)
  if (dist_bin == 0 || dist_bin > length(dist_breaks) - 1) next
  
  colname <- paste0("X", dist_breaks[dist_bin], ".", dist_breaks[dist_bin + 1], "..", det_time)
  out_mat[i, colname] <- 1
}

# Convert to data frame
OR2020_cnt <- as.data.frame(out_mat)

# Optional: retain row IDs from original
# Combine the observation ID as the first column
OR2020_cnt <- cbind(Unique_Observation_ID = OR2020$Unique_Observation_ID, OR2020_cnt)

write.table(OR2020_cnt, file = "/Users/shenf/Documents/Bayes_Class/Final_project/OR2020_cnt.csv", sep = ",", row.names = FALSE, fileEncoding = "UTF-8")
