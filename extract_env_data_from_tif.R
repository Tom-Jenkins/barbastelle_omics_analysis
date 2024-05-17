# ============================ #
#
# Extract Environmental Data from TIF
#
# Data source:
# ./data/CHELSA_climate_data
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(readr)
library(terra)
library(stringr)
library(sf)
library(ggplot2)
library(ggspatial)

# link to raw data
# https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2Ft
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# cmi_mean = Average monthly climate moisture index over 1 year

# Read in coordinates
coords <- read_csv("sample_coordinates.csv")

# Create a SpatVector object from coordinates
coords_sv <- vect(coords, geom = c("LON","LAT"), crs = "epsg:4326")
coords_sv

# Read in 1981-2010 TIFs for each variable
tif_paths <- list.files("./data/CHELSA_climate_data", "1981-2010*", full.names = TRUE)
tif_list <- lapply(tif_paths, rast)
names(tif_list) <- str_extract(tif_paths, "bio\\d+_[0-9][0-9][0-9][0-9]|cmi_[a-z]+_[0-9][0-9][0-9][0-9]")

# Extract data for 1981-2010
tif_data <- lapply(tif_list, function(x) terra::extract(x, coords_sv)[[2]])

# Convert list to data frame
tif_df <- list2DF(tif_data)

# Test multicollinearity
psych::pairs.panels(tif_df, scale = TRUE)

# Only keep key variables that are not correlated (>0.70)
vars_keep <- c("bio5","bio12")
tif_df_noncor <- dplyr::select(tif_df, contains(vars_keep))
head(tif_df_noncor)

# Calculate Variance Inflation Factor
usdm::vif(dplyr::select(tif_df_noncor, contains("1981")))

# Test multicollinearity
psych::pairs.panels(dplyr::select(tif_df_noncor, contains("1981")), scale = TRUE)

# Extract variable values from future data
future_paths <- list.files("./data/CHELSA_climate_data", "2071-2100*", full.names = TRUE)
future_list <- lapply(future_paths, rast)
names(future_list) <- str_extract(future_paths, "bio[0-9]*_\\d{4}-\\d{4}_ukesm1-0-ll_ssp\\d{3}")
future_data <- lapply(future_list, function(x) terra::extract(x, coords_sv)[[2]])
future_df <- list2DF(future_data)

# Add env values to data.frame
coords <- cbind(coords, tif_df_noncor, future_df)
head(coords)

# Export new coordinates file with climate data
write_csv(coords, "sample_coordinates_env.csv")

#--------------#
# Plot maps
#--------------#

# # Study area bounding box for Europe
# europe_bbox <- ext(c(-15, 15, 35, 60))
# 
# # Crop rasters
# bio5_crop <- terra::crop(tif_bio5, europe_bbox)
# plot(bio5_crop)
# bio18_crop <- terra::crop(tif_bio18, europe_bbox)
# plot(bio18_crop)
# 
# # Load coastlines
# library(mapmixture)
# world <- st_read(system.file("extdata", "world.gpkg", package = "mapmixture"))
# europe <- st_crop(world, europe_bbox)
# 
# # Plot
# plot(bio5_crop, axes = "none", legend = "left")
# plot(coords_sv, col = "black", cex = 0.8, alpha = 0.8, add = TRUE)  

