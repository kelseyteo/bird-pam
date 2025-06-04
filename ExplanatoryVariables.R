#################################FOREST COVER##############################
library(raster)
library(terra)
library(sf)
library(landscapemetrics)
library(vegan)
library(ggplot2)
setwd("/Users/kelsey/desktop/")
#Load Map
epping <- rast("testLCM.tif")
print(epping) 
#Upload Coordinates
epping_sites <- read.csv("Epping_Coordinates.csv", header = TRUE)
epping_sites <- st_as_sf(epping_sites, coords=c('Longitude', 'Latitude'), crs=4326)
str(epping_sites)
print(epping_sites)
#clean the data
# Transform the CRS to match the raster (EPSG:27700)
epping_sites <- st_transform(epping_sites, crs = st_crs(epping))
# Get the bounding box and convert it to a polygon
epping_sites_region <- st_as_sfc(st_bbox(epping_sites))

# Buffer the region by 0.1 degrees (approx. 10 km)
#epping_sites_region <- st_buffer(epping_sites_region, 0.1)
epping_sites_region <- st_buffer(epping_sites_region, 10000)

# Crop the raster
epping_sites_landcover <- crop(epping, epping_sites_region)

unique(values(epping))  # Check before cropping
# Extract the first layer (assuming it contains land cover classification)
lcm <- epping[[1]]  
print(lcm)
# Plot to visualize
plot(lcm)
unique(values(lcm))
# Define forest class codes
forest_classes <- c(1, 2)

# Create binary map: 1 = forest, 0 = non-forest
binary_forest <- classify(lcm, cbind(forest_classes, 1), others = 0)

# Plot the binary forest map
plot(binary_forest, col = c("white", "green"), main = "Binary Forest Map")

# Save output
writeRaster(binary_forest, "Binary_Forest_Map.tif", overwrite=TRUE)
epping_sites_forest <- rast("Binary_Forest_Map.tif")
print(epping_sites_forest)
epping_sites_utm23S <- st_transform(epping_sites, 27700)
# This takes a little while to run!
epping_sites_forest_utm23S <- project(epping_sites_forest, "epsg:27700", res=25, method='near') #landcover class is categorical, we HAVE to use method=near, resolution of 30m
plot(epping_sites_forest_utm23S)
plot(st_geometry(epping_sites_utm23S), add=TRUE)

#50M
lsm <- sample_lsm(epping_sites_forest_utm23S, epping_sites_utm23S, 
                  shape = "circle", size = 50, 
                  plot_id = epping_sites_utm23S$Site, 
                  what = c('lsm_c_pland'))
lsm_sil_forest <- subset(lsm, class == 1, select=c(metric, value, plot_id))
# name - the local landscape details are recorded
names(lsm_sil_forest)[2] <- 'C50'
# Reshape that dataset so that each metric has its own column
lsm_sil_forest <- reshape(data.frame(lsm_sil_forest), direction='wide', 
                          timevar='metric', idvar='plot_id')
head(lsm_sil_forest, n=11) #outputs forest cover for each site

#REPEAT UNTIL 800M


######################################NDSI########################################
# Load required packages 
install.packages("tuneR")
install.packages("soundecology")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("seewave")
install.packages("patchwork")
library(tuneR) # For reading audio files
library(soundecology) # For audio analysis functions
library(seewave) # Another package for audio analysis functions
library(dplyr) # To manipulation data
library(ggplot2) # For plotting results
library(tidyr) # For changing the shape of the datasets
library(patchwork) # For combining plots


#Site 4
folder_4 <- "/Users/kelsey/Desktop/EPPING DATA/Birds Site 4" #upload folder of audio files
wav_files_4 <- list.files(path = folder_4, pattern = "\\.(wav|WAV)$", full.names = TRUE)
ndsi_results_4 <- data.frame(File = character(), 
                             NDSI = numeric(),
                             Biophony = numeric(),
                             Anthrophony = numeric(),
                             stringsAsFactors = FALSE)
for (file in wav_files_4) {
  sound <- readWave(file)
  ndsi_result_4 <- ndsi(sound)
  
  ndsi_value_4 <- ndsi_result_4$ndsi_left  # Mono files
  biophony <- ndsi_result_4$biophony_left
  anthrophony <- ndsi_result_4$anthrophony_left
  
  ndsi_results_4 <- rbind(ndsi_results_4, data.frame(File = basename(file), 
                                                     NDSI = ndsi_value_4, 
                                                     Biophony = biophony, 
                                                     Anthrophony = anthrophony))
}
#export to csv
write.csv(ndsi_results_4, "/Users/kelsey/Desktop/EPPING DATA/ndsi_results_4.csv", row.names = FALSE)

#DO THE SAME FOR EVERY SITE