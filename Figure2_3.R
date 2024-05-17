# ============================ #
#
# Barbastelle Figure 2
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ade4)
library(mapmixture)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(patchwork)
library(ggplotify)
library(grid)

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

#--------------#
# Figure 2A: All Genes ####
#--------------#

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
Fig2A <- scatter_plot(
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
  plot_title = "PCA using VST counts of DEGs"
)+
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 14),
  )
Fig2A


#--------------#
# Figure 2A: Heatmap ####
#--------------#

# Set a color palette
heat_cols <- brewer.pal(6, "YlOrRd")

# Create a data frame to annotate heatmap
annotation <- data.frame(
  Temperature = metadata$Temperature,
  Country = metadata$Country,
  row.names = metadata$Sample
)

# Annotation colours list
ann_cols <- list(
  Temperature = c(Low = colours[1], Mid = colours[2], High = colours[3]),
  Country = c(England = "grey90", Portugal = "grey50", Spain = "grey25")
)

# Run pheatmap using metadata for annotation
Fig2A <- pheatmap(
  mat = counts_DEGs,
  annotation_col = annotation,
  annotation_colors = ann_cols,
  scale = "row",
  border_color = "grey30",
  # cellheight = 7,
  # cellwidth = 15,
  cluster_rows = FALSE, 
  show_rownames = FALSE,
  annotation_legend = TRUE,
  labels_col = metadata$site_labels,
  color = heat_cols,
  angle_col = 90,
)

# Remove legend titles
Fig2A$gtable$grobs[[6]]$children[[1]] = textGrob("")
Fig2A$gtable$grobs[[6]]$children[[4]] = textGrob("")


#--------------#
# Figure 2B ####
#--------------#

# Create a new data frame for plotting all DEG normalised counts
counts_DEGs_pt <- counts_DEGs_t

# Add columns for Temperature, Sample and Reproduction
rownames(counts_DEGs_pt) == metadata$Sample
counts_DEGs_pt$Temperature <- metadata$Temperature
counts_DEGs_pt$Sample <- metadata$Sample
counts_DEGs_pt$Reproduction <- metadata$Reproduction

# Pivot data to long format
counts_DEGs_pt <- pivot_longer(
  data = counts_DEGs_pt,
  cols = 1:(ncol(counts_DEGs_pt)-3),
  names_to = "Gene",
  values_to = "Count"
)
counts_DEGs_pt

# Plot using ggplot2
Fig2B <- ggplot(data = counts_DEGs_pt)+
  geom_jitter(aes(x = Gene, y = Count, fill = Temperature), size = 3, shape = 21, width = 0.2)+
  scale_fill_manual(values = colours)+
  # scale_y_log10()+
  ylab("VST counts")+
  xlab("Differentially Expressed Genes")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    # plot.title = element_text(size = 12),
  )+
  guides(fill = guide_legend(override.aes = list(size = 8)))
Fig2B

#--------------#
# Figure 2 Export ####
#--------------#

# Layout design
layout <- "
  A
  B
"

# Plot layout
plt_list = list(
  # Fig2A+ labs(tag = "A")+ theme(legend.position = "none"),
  # as.ggplot(Fig2B)+ labs(tag = "B"),
  # Fig2C+ labs(tag = "C")+ theme(plot.tag = element_text(vjust = -5))
  as.ggplot(Fig2A)+ labs(tag = "A"),
  Fig2B+ labs(tag = "B")
)
Figure2 <- wrap_plots(plt_list, design = layout)
Figure2

# Export
ggsave(plot = Figure2, filename = "Figure2.jpg", width = 15, height = 10, dpi = 600)


#--------------#
# Figure 3 ####
#--------------#

# Function to plot normalised counts
plot_counts <- function(df, gene) {
  
  # subset gene in df
  gene_df <- dplyr::filter(df, Gene == gene)
  
  # plot
  ggplot(data = gene_df, aes(x = Temperature, y = Count))+
    geom_violin(aes(fill = Temperature), scale = "count", alpha = 0.70, colour = NA)+
    scale_fill_manual("Temperature: ", values = c("#91bfdb","#f7f7f7","#fc8d59"))+
    geom_jitter(aes(shape = Reproduction), width = 0.1, fill = "black", colour = "white", size = 3)+
    scale_shape_manual("Reproduction: ", values = c(21,24,22), labels = c("Female (RE)","Female (NR)","Male"))+
    xlab("Temperature")+
    ylab("VST counts")+
    ggtitle(gene)+
    # theme_bw()+
    theme(
      plot.title = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title.y = element_text(size = 15),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.position = "top",
      legend.key = element_rect(fill = NA)
    )+
    guides(fill = guide_legend(override.aes = list(colour = "black")))+
    guides(shape = guide_legend(override.aes = list(size = 5)))
}

# Change labels to Female (L/PL) Female (NR) Male 

# Custom
# plot_counts(counts_DEGs_pt, "IP6K1")
# ggsave("Figure S1 IP6K1.jpg", width = 10, height = 5, dpi = 600)

# DNAJC4 
(DNAJC4 <- plot_counts(counts_DEGs_pt, "DNAJC4"))

# UCP2 
(UCP2 <- plot_counts(counts_DEGs_pt, "UCP2"))

# DMTN
(DMTN <- plot_counts(counts_DEGs_pt, "DMTN"))

# RAPGEF5
(RAPGEF5 <- plot_counts(counts_DEGs_pt, "RAPGEF5"))

# AGRN 
(AGRN <- plot_counts(counts_DEGs_pt, "AGRN"))

# MED25 
(MED25 <- plot_counts(counts_DEGs_pt, "MED25"))

# BAG6 
(BAG6 <- plot_counts(counts_DEGs_pt, "BAG6"))

# LRP10 
(LRP10 <- plot_counts(counts_DEGs_pt, "LRP10"))

# MFSD12 
(MFSD12 <- plot_counts(counts_DEGs_pt, "MFSD12"))

# New ggplot2 theme
gg_theme <- theme(
  legend.position = "top",
  legend.title = element_text(size = 20, face = "bold"),
  legend.text = element_text(size = 20),
  plot.title = element_text(size = 20),
)

# Plot layout
plt_list = list(
  UCP2+ gg_theme,
  DNAJC4+ gg_theme+ theme(axis.title.y = element_blank()),
  BAG6+ gg_theme+ theme(axis.title.y = element_blank()),
  DMTN+ gg_theme,
  MED25+ gg_theme+ theme(axis.title.y = element_blank()),
  AGRN+ gg_theme+ theme(axis.title.y = element_blank())
)
Figure3 <- wrap_plots(plt_list, ncol = 3)+
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
Figure3

# Export
ggsave(plot = Figure3, filename = "Figure3.jpg", width = 15, height = 8, dpi = 600)


