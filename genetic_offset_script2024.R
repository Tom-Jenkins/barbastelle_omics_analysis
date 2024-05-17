# ============================ #
#
# Genetic Offset Analysis
#
# Data source:
# ./data/barb_m50g10maf03_2024.lfmm_imputed.lfmm
# ./sample_coordinates_env.csv
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(LEA)
library(readr)
library(data.table)
library(ggplot2)
library(ggtext)
library(dplyr)
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(mapmixture)
library(ade4)

# link to climate data
# https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2Ft
# BIO5 = Max Temperature of Warmest Month
# BIO12 = Annual Precipitation

# Read in genotypes in LFMM format
geno_outlier <- data.table::fread("./data/barb_m3_g90_maf03_feb2024.imputed.lfmm.outliers")

# Read in sample data and convert env data to a matrix
sample_data <- read_csv("sample_coordinates_env.csv")
env_1981_2010 <- sample_data[ , c("bio5_1981","bio12_1981")]
env_2071_2100_ssp126 <- sample_data[ , c("bio5_2071-2100_ukesm1-0-ll_ssp126","bio12_2071-2100_ukesm1-0-ll_ssp126")]
env_2071_2100_ssp370 <- sample_data[ , c("bio5_2071-2100_ukesm1-0-ll_ssp370","bio12_2071-2100_ukesm1-0-ll_ssp370")]
env_2071_2100_ssp585 <- sample_data[ , c("bio5_2071-2100_ukesm1-0-ll_ssp585","bio12_2071-2100_ukesm1-0-ll_ssp585")]

#--------------#
# Plot Geno-Env Correlation ####
#--------------#

# Run PCA
pca_geno <- dudi.pca(geno_outlier, scale = TRUE, scannf = FALSE, nf = 3)

# Correlation between sample PC1,PC2 and sample env data
plot(pca_geno$li$Axis1, env_1981_2010$bio5_1981); cor(pca_geno$li$Axis1, env_1981_2010$bio5_1981)
plot(pca_geno$li$Axis1, env_1981_2010$bio12_1981); cor(pca_geno$li$Axis1, env_1981_2010$bio12_1981)
plot(pca_geno$li$Axis2, env_1981_2010$bio5_1981); cor(pca_geno$li$Axis2, env_1981_2010$bio5_1981)
plot(pca_geno$li$Axis2, env_1981_2010$bio12_1981); cor(pca_geno$li$Axis2, env_1981_2010$bio12_1981)

#--------------#
# Genomic Offset Analysis ####
#--------------#

# Run genetic offset function
genomic_offset_ssp585 <- genetic.offset(
  input = geno_outlier, 
  env = env_1981_2010,
  pred.env = env_2071_2100_ssp585,
  K = 3
)

# Return geometric genomic offset (genetic gap) for each sample location
goffset_ssp585 <- round(genomic_offset_ssp585$offset, digit = 2) 

# Return RONA for each sample location
rona_ssp585 <- round(genomic_offset_ssp585$distance, digit = 2) 

# Add both to sample data frame
sample_data$offset_ssp585 <- goffset_ssp585
# sample_data$rona_ssp585 <- rona_ssp585

# Plot RONA vs. the genetic gap 
plot(genomic_offset_ssp585$offset, genomic_offset_ssp585$distance)

# Re-run but with scaled variables
scaled <- genetic.offset(
  input = geno_outlier, 
  env = env_1981_2010,
  pred.env = env_2071_2100_ssp585,
  K = 3,
  scale = TRUE
)

# Scaling does not change genetic offsets
plot(genomic_offset_ssp585$offset, scaled$offset)                           

# But scaling is useful for evaluating the relative importance of environmental variables
# Two dimensions in environmental space have influence on the genetic offset
barplot(scaled$eigenvalues, col = "orange", xlab = "Axes", ylab = "Eigenvalues")

# The loadings for the first two variables indicate their relative contribution to local adaptation
scaled$vectors[,1:2]

# Run genetic offset function
genomic_offset_ssp126 <- genetic.offset(
  input = geno_outlier, 
  env = env_1981_2010,
  pred.env = env_2071_2100_ssp126,
  K = 3
)
goffset_ssp126 <- round(genomic_offset_ssp126$offset, digit = 2) 
sample_data$offset_ssp126 <- goffset_ssp126

#--------------#
# Basemap ####
#--------------#

# CRS
CRS = 3035

# Europe basemap
basemap <- ne_countries(scale = "large", continent = "Europe")[, c("admin")]
basemap <- st_transform(basemap, crs = CRS)

# Boundary
bbox <- mapmixture::transform_bbox(c(xmin = -12, xmax = 8, ymax = 38, ymin = 54), CRS)

# Convert sample data frame to sf object
sample_sf <- st_transform(st_as_sf(sample_data, coords = c("LON", "LAT"), crs = 4326), crs = CRS)

# Add changes in climate variables
sample_sf$bio5_ssp126_change <- sample_sf$`bio5_2071-2100_ukesm1-0-ll_ssp126` - sample_sf$bio5_1981
sample_sf$bio12_ssp126_change <- sample_sf$`bio12_2071-2100_ukesm1-0-ll_ssp126` - sample_sf$bio12_1981
sample_sf$bio5_ssp585_change <- sample_sf$`bio5_2071-2100_ukesm1-0-ll_ssp585` - sample_sf$bio5_1981
sample_sf$bio12_ssp585_change <- sample_sf$`bio12_2071-2100_ukesm1-0-ll_ssp585` - sample_sf$bio12_1981

# Summarise data by location and take the centroid of the lat and lon for each geometry
location_sf <- sample_sf |>
  mutate(location = word(sample_data$LOCATION, start = 1, sep = ",")) |>
  group_by(location) |>
  dplyr::summarise(
    offset_ssp126_mean = mean(offset_ssp126),
    offset_ssp585_mean = mean(offset_ssp585),
    bio5_ssp585_change_mean = mean(bio5_ssp585_change),
    bio5_ssp585_change_sd = sd(bio5_ssp585_change),
    bio5_ssp126_change_mean = mean(bio5_ssp126_change),
    bio5_ssp126_change_sd = sd(bio5_ssp126_change),
    bio12_ssp585_change_mean = mean(bio12_ssp585_change),
    bio12_ssp585_change_sd = sd(bio12_ssp585_change),
    bio12_ssp126_change_mean = mean(bio12_ssp126_change),
    bio12_ssp126_change_sd = sd(bio12_ssp126_change)
  ) |>
  st_centroid()

# Theme
gg_theme <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "#f0f0f0"),
  plot.title = element_text(size = 10)
)

#--------------#
# Figure: Predicted Climate Changes ####
#--------------#

# BIO5 plot ssp126
temp1 <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = location_sf, aes(fill = bio5_ssp126_change_mean), size = 5, shape = 21, colour = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  scale_fill_gradient(name = "Temperature<br>change (<sup>o</sup>C)", 
                      low = "#ffeda0", high = "#f03b20",
                      limits = c(3.0,10.5), breaks = seq(4,10,2))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Max. Temperature of Warmest Month – SSP126")+
  gg_theme+
  theme(legend.title = ggtext::element_markdown())
temp1

# BIO5 plot ssp585
temp2 <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = location_sf, aes(fill = bio5_ssp585_change_mean), size = 5, shape = 21, colour = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  scale_fill_gradient(name = "Temperature<br>change (<sup>o</sup>C)", 
                      low = "#ffeda0", high = "#f03b20",
                      limits = c(3.0,10.5), breaks = seq(4,10,2))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Max. Temperature of Warmest Month – SSP585")+
  gg_theme+
  theme(legend.title = ggtext::element_markdown())
temp2

# BIO12 plot ssp126
rain1 <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = location_sf, aes(fill = bio12_ssp126_change_mean), size = 5, shape = 21, colour = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  scale_fill_gradient2(name = "Precipitation\nchange (mm)",
                       low = "#8c510a", mid = "white", high = "#01665e",
                       midpoint = 0, limits = c(-166,73), breaks = seq(-150,50,50))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Total Annual Precipitation – SSP126")+
  gg_theme
rain1

# BIO12 plot ssp585
rain2 <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = location_sf, aes(fill = bio12_ssp585_change_mean), size = 5, shape = 21, colour = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  scale_fill_gradient2(name = "Precipitation\nchange (mm)", low = "#8c510a", mid = "white", high = "#01665e",
                       midpoint = 0, limits = c(-166,73), breaks = seq(-150,50,50))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Total Annual Precipitation – SSP585")+
  gg_theme
rain2

#--------------#
# Figure: Genomic Offset ####
#--------------#

# Genomic offset map ssp126
offset1 <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = location_sf, aes(fill = offset_ssp126_mean), size = 5, shape = 21, colour = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  # scale_fill_gradient(name = "Genomic Offset", low = "green3", high = "deeppink",
  #                     limits = c(0,0.45), breaks = seq(0,0.40,0.10))+
  scale_fill_viridis_c(name = "Genomic Offset", option = "viridis",
                       limits = c(0,0.45), breaks = seq(0,0.40,0.10))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Geometric Genomic Offset – SSP126")+
  gg_theme
offset1

# Genomic offset map ssp585
offset2 <- ggplot()+
  geom_sf(data = basemap)+
  geom_sf(data = location_sf, aes(fill = offset_ssp585_mean), size = 5, shape = 21, colour = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  # scale_fill_gradient(name = "Genomic Offset", low = "green3", high = "deeppink",
  #                     limits = c(0,0.45), breaks = seq(0,0.40,0.10))+
  scale_fill_viridis_c(name = "Genomic Offset", option = "viridis",
                       limits = c(0,0.45), breaks = seq(0,0.40,0.10))+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Geometric Genomic Offset – SSP585")+
  gg_theme
offset2


#--------------#
# Figure: Composer ####
#--------------#

# Load patchwork
library(patchwork)

# Layout design
plt_list = list(
  temp1 + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(tag = "A"),
  temp2 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  rain1 + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(tag = "B"),
  rain2 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  offset1 + theme(legend.position = "none") + labs(tag = "C"),
  offset2 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
)
fig <- wrap_plots(plt_list, nrow = 3, ncol = 2)

# Export figure
ggsave(plot = fig, filename = "../Figure5.jpg", width = 10, height = 13, units = "in", dpi = 600)

