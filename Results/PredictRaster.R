#---
#title: "Integration Scrub Jay"
#author: "Betty Shen"
#date: "2025-05-24"
#output: html_document
#---
#Goal: the goal of the script is to predict density of Scrub Jay in 2019 using IDS model

#Load library
#install.packages("unmarked")
library(unmarked)
#packageVersion("unmarked")
library(dplyr)
library(raster)
library(sf)
library(terra)


# Load raster

raster <- readRDS("/nfs/stak/users/shenf/hpc-share/Test/xy2019_Landsat_PRISM_EVI.rds")

print(names(raster))
# Rename columns
raster <- raster %>%
  rename(
    SR_B1_median = Landsat_PRISM_EVI2019_clipped_1,
    SR_B5_median = Landsat_PRISM_EVI2019_clipped_2,
    EVI_median = Landsat_PRISM_EVI2019_clipped_3,
    ppt_median = Landsat_PRISM_EVI2019_clipped_4,
    tmean_median = Landsat_PRISM_EVI2019_clipped_5
  )
# Only keep the columns that match the model's covariates
newdata <- raster[,c(3:7)]
names(newdata)

# Scale the predictors
newdata_scaled <- newdata
newdata_scaled$EVI_median_scaled <- scale(newdata$EVI_median)
newdata_scaled$SR_B1_median_scaled <- scale(newdata$SR_B1_median)
newdata_scaled$SR_B5_median_scaled <- scale(newdata$SR_B5_median)
newdata_scaled$ppt_median_scaled <- scale(newdata$ppt_median)
newdata_scaled$tmean_median_scaled <- scale(newdata$tmean_median)

# Use the scaled data for prediction
newdata <- newdata_scaled[, c("EVI_median_scaled", "SR_B1_median_scaled", "SR_B5_median_scaled", 
                             "ppt_median_scaled", "tmean_median_scaled")]
names(newdata)

print("Load Scrub Jay IDS model")
# Load IDS model
IDS.model <- readRDS("/nfs/stak/users/shenf/hpc-share/IntegrationModel/TrainedModel/SpeciesSpecif/California Scrub-Jay.rds")

# Predict density per pixel
print("Predict density per pixel")
# Predict values from the IDS model in chunks of 1 million rows
chunk_size <- 1000000
n_rows <- nrow(newdata)
n_chunks <- ceiling(n_rows / chunk_size)

# Initialize empty list to store predictions
all_preds <- list()

# Loop through chunks
for(i in 1:n_chunks) {
  start_row <- ((i-1) * chunk_size) + 1
  end_row <- min(i * chunk_size, n_rows)
  
  # Get chunk of data
  chunk_data <- newdata[start_row:end_row,]
  
  # Predict for this chunk
  chunk_preds <- predict(IDS.model, newdata = chunk_data, type = "lam")
  
  # Store predictions
  all_preds[[i]] <- chunk_preds
  
  # Print progress
  print(paste("Completed chunk", i, "of", n_chunks))
}

# Combine all predictions
preds <- do.call(rbind, all_preds)

# Bind predictions back to xy
#raster$prediction <- preds$Predicted

# Save predictions
saveRDS(preds, "/nfs/stak/users/shenf/hpc-share/Test/FullRaster_ScrubJay_preds.rds")

preds <- readRDS("/nfs/stak/users/shenf/hpc-share/Test/FullRaster_ScrubJay_preds.rds")
str(preds)
range(preds$Predicted)
# We need to downscale the predicted values from 1 km^2 to 90 m^2 
preds_scale <- preds  # Create a copy of the original dataframe
str(preds_scale)
preds_scale[c("Predicted", "lower", "upper")] <- preds[c("Predicted", "lower", "upper")] * (9/10000)

# Merge with original raster file (update to use preds_scale)
raster <- cbind(raster, preds_scale)
str(raster)

# Convert to a new raster
pred_rast <- rast(raster[,c("x","y","Predicted")], type = "xyz")
saveRDS(pred_rast, "/nfs/stak/users/shenf/hpc-share/Test/pred_rast_ScrubJay2019.rds")


# Read in the original raster
hist(values(pred_rast), breaks=100, main="Distribution of Predicted Values")
range(pred_rast)

cols <- hcl.colors(100, "Viridis")

# Cap values above 0.7 (after downscaling)
pred_rast_capped <- pred_rast
pred_rast_capped[pred_rast_capped > 0.3] <- 0.3 

png("/nfs/stak/users/shenf/hpc-share/Test/Raster_ScrubJay2019_scaled.png", width=1200, height=900)
plot(pred_rast_capped, zlim=c(0, 0.3), col=cols)
print("Plot done")
dev.off()
