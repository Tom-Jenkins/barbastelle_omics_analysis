# ============================ #
#
# Barbastelle Figure 1
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(jpeg)
library(grid)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(dplyr)
library(stringr)
library(mapmixture)
library(readr)
# library(ggrepel)
library(terra)
library(ggspatial)
library(patchwork)
# library(ggnewscale)
library(ade4)

#--------------#
# Figure 1A ####
#--------------#

# CRS
CRS = 3035

# Europe basemap
basemap <- ne_countries(scale = "large", continent = "Europe")[, c("admin")]
basemap <- st_transform(basemap, crs = CRS)

# Boundary
bbox <- mapmixture::transform_bbox(c(xmin = -12, xmax = 8, ymax = 37, ymin = 57), CRS)

# Read in sample metadata csv as sf object
metadata <- read_csv("Table S1_ Samples_metadata.csv") |>
  st_as_sf(x = _, coords = c("Lon","Lat"), crs = 4326) |>
  st_transform(x = _, crs = CRS)
metadata

# Create a new data frame of unique locations
metadata_plt <- dplyr::select(metadata, Location) |> distinct()

# Create new column with only county name
metadata_plt$county <- word(metadata_plt$Location, start = 1, sep = ",")
metadata_plt

# Create new data frame with only unique county (the centroid of each county)
metadata_county <- metadata_plt |>
  group_by(county) |>
  summarise() |>
  st_centroid()

# Add map labels to data frame
metadata_county$map_labels <- substr(metadata_county$county, start = 1, stop = 3)
metadata_county$map_labels <- str_replace(metadata_county$map_labels, "La ", "Rio")
metadata_county

# Column for temperature group
metadata_county$temp <- c(
  "Bed"="Mid",
  "Dar"="RAD Only",
  "Gal"="Low",
  "Jae"="Mid High",
  "Rio"="Low Mid",
  "Not"="Mid",
  "Sab"="High",
  "Sor"="Low",
  "Sus"="Mid",
  "Ter"="Low High",
  "War"="Mid"
)
metadata_county

# Read in land cover raster
# elev_rast <- terra::rast("NE1_HR_LC_SR/NE1_HR_LC_SR.tif")
# elev_crop <- terra::crop(elev_rast, y = c(xmin = -30, xmax = 10, ymax = 30, ymin = 70))
# terra::writeRaster(elev_crop, filename = "NE1_HR_LC_SR/NE1_HR_LC_SR_crop.tif", overwrite = TRUE)
elev_rast <- rast("NE1_50M_SR_W/NE1_50M_SR_W_crop.tif") |> terra::project(x = _, y = str_c("epsg:", CRS))
plot(elev_rast)

# Convert all values below sea level (less than zero) to NA
# elev_rast <- terra::clamp(elev_rast, lower = 0, values = TRUE)

# Extract elevation for all coordinates
# terra::extract(elev_rast, y = metadata_plt) |> mutate(location = metadata_plt$Location)

# Extract XY coordinates
coordinates <- as.data.frame(st_coordinates(metadata_county))
coordinates$site <- metadata_county$map_labels

# Theme
gg_theme <- theme(
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "#f0f0f0"),
  plot.title = element_text(size = 10)
)

# Plot sampling map
Fig1A <- ggplot()+
  layer_spatial(data = elev_rast)+
  scale_fill_gradientn(
    name = "Elevation (m)",
    colours = c("#668E47", "#DAB263", "#A98712", "#865B11", "#84473F"),
    na.value = NA
  )+
  annotation_north_arrow(
    data = basemap, location = "tl", which_north = "true",
    height = unit(0.6, "cm"), width = unit(0.6, "cm"),
    pad_y = unit(0.8, "cm"),
    style = north_arrow_orienteering(text_size = 5)
  )+
  annotation_scale(
    data = basemap, location = "tl",
    width_hint = 0.2, bar_cols = c("black","white"),
    # height = unit(0.5, "cm"),
    text_cex = 0.8
  )+
  annotate("label", x = coordinates[1,]$X+100000, y = coordinates[1,]$Y, label = "Bed", fill = "#f7f7f7", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[2,]$X-100000, y = coordinates[2,]$Y, label = "Dar", fill = "grey", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[3,]$X-100000, y = coordinates[3,]$Y, label = "Gal", fill = "#67a9cf", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[4,]$X-130000, y = coordinates[4,]$Y, label = "Jae16", fill = "#f7f7f7", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[4,]$X+130000, y = coordinates[4,]$Y, label = "Jae18", fill = "#ef8a62", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[5,]$X-150000, y = coordinates[5,]$Y, label = "Rio16", fill = "#f7f7f7", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[5,]$X, y = coordinates[5,]$Y+80000, label = "Rio18", fill = "#67a9cf", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[6,]$X, y = coordinates[6,]$Y+80000, label = "Not", fill = "#f7f7f7", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[7,]$X-100000, y = coordinates[7,]$Y, label = "Sab", fill = "#ef8a62", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[8,]$X+100000, y = coordinates[8,]$Y, label = "Sor", fill = "#67a9cf", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[9,]$X+100000, y = coordinates[9,]$Y, label = "Sus", fill = "#f7f7f7", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[10,]$X+130000, y = coordinates[10,]$Y, label = "Ter16", fill = "#67a9cf", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[10,]$X-130000, y = coordinates[10,]$Y, label = "Ter18", fill = "#ef8a62", size = 3.5, label.padding = unit(0.15, "lines"))+
  annotate("label", x = coordinates[11,]$X-100000, y = coordinates[11,]$Y, label = "War", fill = "#f7f7f7", size = 3.5, label.padding = unit(0.15, "lines"))+
  geom_sf(data = metadata_county, size = 2.5, shape = 21, colour = "white", fill = "black")+
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))+
  xlab("Longitude")+
  ylab("Latitude")+
  gg_theme
Fig1A

#--------------#
# Figure 1B ####
#--------------#

# Import jpeg and convert to ggplot object
barb_jpeg <- readJPEG("barbastelle.jpg") |> rasterGrob(image = _, interpolate = TRUE)

# Plot image
Fig1B <- ggplot()+
  annotation_custom(barb_jpeg, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_void()
Fig1B

#--------------#
# Figure 2A: All Genes ####
#--------------#

# Read in sample metadata csv and subset RNA samples
metadata <- read_csv("Table S1_ Samples_metadata.csv") |>
  dplyr::filter(RNA == "Y") |>
  arrange(Sample)
metadata

# Create column for site labels
metadata$site_labels <- str_extract(metadata$Location, "^.{3}")
metadata$site_labels <- str_replace(metadata$site_labels, "La ", "Rio")

# Assign a condition grouping to measured temperature (celsius)
metadata <- metadata |>
  mutate(.data = _, Temperature = case_when(
    Site_Temp_C >= 9.9 & Site_Temp_C <= 11.3 ~ "Low",
    Site_Temp_C >= 13 & Site_Temp_C <= 14.5 ~ "Mid",
    Site_Temp_C >= 15 & Site_Temp_C <= 16.4 ~ "High",
  )) |>
  mutate(.data = _, Temperature = factor(Temperature, levels = c("Low","Mid","High")))

# Temperature groups
dplyr::count(metadata, Temperature)

# Read in PCA results
# counts <- read.csv("Gene_expression_analysis/data/normalised_counts.csv", row.names = 1)
counts <- read.csv("Gene_expression_analysis/data/vst_counts.csv", row.names = 1)

# Colour palette
colours <- c("#67a9cf","#f7f7f7","#ef8a62")

# Extract all genes and transpose data frame
counts_t <- as.data.frame(t(counts))

# Check sample order
rownames(counts_t) == metadata$Sample

# Re-run PCA
pca1 <- dudi.pca(counts_t, scale = TRUE, scannf = FALSE, nf = 10)

# Analyse how much percent of genetic variance is explained by each axis
percent <- round(pca1$eig/sum(pca1$eig)*100, digits = 1)

# Plot PCA results and colour by temperature
scatter_plot(
  dataframe = as.data.frame(pca1$li),
  type = "labels",
  group_ids = metadata$Temperature,
  labels = metadata$site_labels,
  label.padding = unit(0.20, "lines"),
  size = 3,
  axes = c(1,2),
  percent = percent,
  colours = colours,
  plot_title = "PCA with all normalised gene counts"
)

#--------------#
# Figure 2A: DEGs ####
#--------------#

# Read in DEGs
DEGs <- read.csv("Gene_expression_analysis/DEGs.csv") |> distinct()

# Add a nonDEG to DEGs (only for visualation)
# nonDEGs <- "IP6K1"
# DEGs <- rbind(DEGs, nonDEGs)

# Extract only counts for DEGs
counts_DEGs <- counts[which(rownames(counts) %in% DEGs$DEG), ]

# Transpose data frame and add gene names to columns
counts_DEGs_t <- as.data.frame(t(counts_DEGs))

# Check sample order
rownames(counts_DEGs_t) == metadata$Sample

# Re-run PCA
pca2 <- dudi.pca(counts_DEGs_t, scale = TRUE, scannf = FALSE, nf = 10)

# Analyse how much percent of genetic variance is explained by each axis
percent_DEGs <- round(pca2$eig/sum(pca2$eig)*100, digits = 1)

# Plot PCA results and colour by temperature
Fig1C <- scatter_plot(
  dataframe = as.data.frame(pca2$li),
  type = "labels",
  group_ids = metadata$Temperature,
  labels = metadata$site_labels,
  label.padding = unit(0.20, "lines"),
  size = 3,
  axes = c(1,2),
  xlab = "PC",
  ylab = "PC",
  percent = percent_DEGs,
  colours = colours,
  plot_title = "PCA of DEGs identified in temperature model",
)+
  scale_y_continuous(position = "right")+
  scale_fill_manual(
    name = "Temperature",
    labels = c("Low", "Mid", "High"),
    values = colours
  )+
  theme(
    # plot.title = element_text(size = 15),
    plot.title = element_blank(),
    axis.title = element_text(size = 13),
    legend.position = "left",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(2, "lines"),
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 4, label = "  ")))
Fig1C


#--------------#
# Figure 1 Export ####
#--------------#

# Layout design
layout <- "
  AABB
  AACC
  AACC
"

# Plot layout
plt_list = list(
  Fig1A+ labs(tag = "A"),
  Fig1B+ labs(tag = "B"),
  Fig1C+ labs(tag = "C")
)
Figure1 <- wrap_plots(plt_list, design = layout)

# Export
ggsave(plot = Figure1, filename = "Figure1.jpg", width = 13, height = 7, dpi = 600)


#--------------#
# Figure S1 ####
#--------------#

# Import RData file containing PCA results
load("Genotype_environment_association_analysis/pca_popgen_46230snps.RData")

# Analyse how much percent of genetic variance is explained by each axis
percent <- pca_geno$eig/sum(pca_geno$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)")

# Read in sample metadata and subset RAD samples
metadata_rad <- read_csv("Table S1_ Samples_metadata.csv") |>
  dplyr::filter(!is.na(ddRAD)) |>
  mutate(Order = as.numeric(str_remove(ddRAD, "Y"))) |>
  dplyr::arrange(Order)
metadata_rad

# Dartmoor point
metadata_dart <- pca_geno$li[which(metadata_rad$Location == "Dartmoor, Devon"),]

# Character vector of site IDs
site_ids <- factor(c(rep("Sussex", 23), "England", rep("Sussex", 12), rep("England", 20), rep("Iberia", 38), "England"), levels = c("England","Sussex","Iberia"))

# Character vector of labels
labels <- str_extract(metadata_rad$Location, "^.{2}")

# Plot PCA
FigS1 <- scatter_plot(
  dataframe = pca_geno$li, group_ids = site_ids, type = "labels", labels = labels,
  percent = round(percent, digits = 1),
  xlab = "PC",
  ylab = "PC",
  colours = c("#998ec3","#d8daeb","#f1a340"), 
  size = 3, alpha = 0.95, label.padding = unit(0.15, "lines"),
)+
  annotate(geom = "label", x = metadata_dart$Axis1, metadata_dart$Axis2, label = "Da",
           fill = "#998ec3", size = 3, label.padding = unit(0.15, "lines"))+
  annotate(geom = "text", x = 57, y = 100, label = "PCA using 46,230 SNPs", size = 3)+
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.35, "lines"),
    legend.spacing.x = unit(0.25, "cm"),
    axis.title = element_text(size = 10)
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 4, label = "  ")))
FigS1
ggsave(plot = FigS1, filename = "Figure S1_Population_Structure.jpg", width = 10, height = 7, units = "in")
