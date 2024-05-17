# ============================ #
#
# Barbastelle Differential Gene Expression Analysis
#
# Tissue: Blood
# Data: quant.sf outputs from `salmon quant <args>` command
# Date: Feb 2024
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(tximport)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(DESeq2)
library(DEGreport)
library(apeglm)
library(ashr)
library(pheatmap)
library(rlang)
library(ggrepel)

# Source utility functions
source("utils.R")

# ----------------- #
# Salmon Quant results
# ----------------- #

# Read in gene IDs for each transcript
gene_ids <- read_tsv(
  file = "data/id_transcript_gene.tsv",
  show_col_types = FALSE,
  col_names = c("transcript","gene")
)

# Salmon files to analyse
files <- list.files("./data/Quant_fastp/", pattern = "*.sf", recursive = T, full.names = T)
files <- str_subset(files, "_BL")

# Extract sample names from files
sample_names = str_extract(files, "B[A-Za-z0-9]*_[A-Z][A-Z]")

# Read in count data using tximport
txi <- tximport(
  files = files,
  type = "salmon",
  tx2gene = gene_ids,
  countsFromAbundance = "lengthScaledTPM"
)

# ----------------- #
# Sample metadata #### 
# ----------------- #

# Read in sample metadata
metadata <- read_csv("../Samples_metadata.csv", show_col_types = FALSE) |>
  dplyr::filter(RNA == "Y") |>
  arrange(.data = _, Sample) |>
  dplyr::select(.data = _, Sample, Location, Country, Sex, Reproduction, Site_Temp_C)
metadata

# Check order of file import names matches order of metadata names
metadata$Sample == sample_names

# Reproduction groups
dplyr::count(metadata, Reproduction)

# Assign a condition grouping to measured temperature (celsius)
sort(metadata$Site_Temp_C)
metadata <- metadata |>
  mutate(.data = _, Temperature = case_when(
    Site_Temp_C >= 9.9 & Site_Temp_C <= 11.3 ~ "low",
    Site_Temp_C >= 13 & Site_Temp_C <= 14.5 ~ "mid",
    Site_Temp_C >= 15 & Site_Temp_C <= 16.4 ~ "high",
  ))

# Temperature groups
dplyr::count(metadata, Temperature)

# Prepare coldata matrix input to DeSeq2
coldata_mat <- metadata |>
  dplyr::select(Temperature, Country, Sex, Reproduction) |>
  as.matrix(x = _)
rownames(coldata_mat) <- metadata$Sample
head(coldata_mat)

# ----------------- #
# Plot count data
# ----------------- #

# Subset count data and round to nearest integer
# count_data <- txi$counts |>
#   round(x = _, digits = 0) |>
#   as.data.frame(x = _) |>
#   set_names(sample_names)

# Plot counts for a single sample
# sample_to_plt <- sample_names[1]
# ggplot(count_data)+
#   geom_histogram(aes(x = !!as.name(sample_to_plt)), stat = "bin", bins = 200)+
#   xlab("Raw expression counts")+
#   ylab("Number of genes")+
#   scale_x_continuous(labels = scales::comma)

# Plot mean versus variance for Curro do Pereiro, Galicia blood samples
# sample_sub <- c("B170710Bba1_BL","B170710Bba2_BL","B170713Bba1_BL")
# Bba_mean <- apply(count_data[,sample_sub], 1, mean)
# Bba_variance <- apply(count_data[,sample_sub], 1, var)
# ggplot(data = data.frame(Bba_mean, Bba_variance))+
#   geom_point(aes(x=Bba_mean, y=Bba_variance))+
#   scale_y_log10(limits = c(1,1e10))+
#   scale_x_log10(limits = c(1,1e10))+
#   geom_abline(intercept = 0, slope = 1, color="red")+
#   xlab("Mean")+
#   ylab("Variance")

# Variance across biological replicates greater than the mean
# Indicates that our data do not fit a Poisson distribution (mean == variance)
# Therefore the negative binomial model is a better fit (mean < variance)

# ----------------- #
# Model parameters
# ----------------- #

# Set p-value threshold
alpha <- 0.01

# Log2 fold change (LFC) threshold (1.5 = 0.58 effect size)
# lfc_threshold <- log2(1.5)
lfc_threshold <- 1.5

# ----------------- #
# Reproduction model ####
# ----------------- #

# Construct a DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi, coldata_mat, ~ Reproduction)

# Reorder factor levels
dds$Reproduction <- factor(dds$Reproduction, levels = c("Male", "Female_nonreproductive", "Female_lac_postlac"))
dds$Reproduction

# Estimate range of size factors using Median of ratios normalisation method
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
dds <- estimateSizeFactors(dds)
sort(sizeFactors(dds))
# plot(sizeFactors(dds))

# # Remove outlier samples
# deceased_bat <- c("B680_BL")
# lowSizeFactorRatios <- c("B170710Bba1_BL","B582_BL","B158_BL","B258_BL")
# highSizeFactorRatios <- c("B296_BL","B691_BL","B170713Bba1_BL","B234_BL","B294_BL","B254_BL")
# outlier_samples <- c(deceased_bat, lowSizeFactorRatios, highSizeFactorRatios)
# dds <- dds[, !(colnames(dds) %in% outlier_samples)]
# dds
# 
# # Estimate size factors
# dds <- estimateSizeFactors(dds)
# sort(sizeFactors(dds))


# Print number of samples per group
table(dds[["Reproduction"]])

# QC: Gene-level
# Pre-filter: omit genes with zero counts in all samples
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]

# QC: Gene-level
# Pre-filter: omit genes that have less than 10 counts for the smallest group size
# E.g. if minGroupSize=5, five or more samples must have at least 10 counts
minGroupSize <- min(table(dds[["Reproduction"]]))
keep <- rowSums(counts(dds) >= 10) >= minGroupSize
dds <- dds[keep, ]
dds

# Total number of counts per sample
# colSums(counts(dds))

# Total number of normalized counts per sample
# colSums(counts(dds, normalized = TRUE))

# Extract normalised counts
normalised_counts <- counts(dds, normalized = TRUE)

# Export normalised counts
# write.csv(normalised_counts, file = "data/normalised_counts.csv", row.names = T)

# Variance Stabilising Transformation (VST)
vst <- varianceStabilizingTransformation(dds, blind = TRUE)
# assay(vst)[1:5, 1:5]

# QC: Sample-level
# Principal components analysis
# pltPCA(vst, metadata, group = "Temperature", axes = c("PC1","PC2"), size = 5)
# pltPCA(vst, metadata, group = "Sex", axes = c("PC1","PC2"), size = 5)
# pltPCA(vst, metadata, group = "Reproduction", axes = c("PC1","PC2"), size = 5)
# pltPCA(vst, metadata, group = "Country", axes = c("PC1","PC2"), size = 5)

# Run analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Coefficients estimated
resultsNames(dds)

# Groups to compare
comparison <- c("Reproduction", "Female_lac_postlac", "Male")
# comparison <- c("Reproduction", "Female_lac_postlac", "Female_nonreproductive")
# comparison <- c("Reproduction", "Female_nonreproductive", "Male")

# Results table
res <- results(dds, contrast = comparison, alpha = alpha)
summary(res)

# Shrink log2 fold changes using adaptive shrinkage estimator (ashr) method
res_shrunken <- lfcShrink(dds, contrast = comparison, res = res, type = "ashr")
summary(res_shrunken)

# Plot normalized counts mean versus the log2 foldchanges for all genes tested
# plotMA(res, alpha = alpha, ylim = c(-5,5), colSig = "red", main = "Results")
# plotMA(res_shrunken, alpha = alpha, ylim = c(-5,5), colSig = "red", main = "Results Shrunken")

# Plot log2 fold change versus adjusted p-values
volcano_plt <- as_tibble(res_shrunken) |>
  mutate(.data = _, gene = rownames(res_shrunken)) |>
  mutate(.data = _, DEG = padj < alpha) |>
  drop_na(data = _) |>
  mutate(.data = _, gene_labels = ifelse(DEG == T, gene, "")) |>
  ggplot(data = _)+
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = DEG), size = 3)+
  xlab("log2 fold change")+
  ylab("-log10 adjusted p-value")+
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = gene_labels))
# volcano_plt

# DEGs table
DEGs <- as_tibble(res_shrunken) |>
  # Add column for gene names
  mutate(.data = _, gene = rownames(res_shrunken)) |>
  # Reorder data frame
  select(.data = _, gene, everything()) |>
  # Remove columns with NAs
  drop_na(data = _) |>
  # Filter genes with adjusted p-values < alpha
  subset(x = _, padj < alpha) |>
  # Filter genes with absolute values > lfc_threshold
  subset(x = _, abs(log2FoldChange) > lfc_threshold) |>
  # Order by adjusted p-value
  arrange(.data = _, padj)
DEGs
# write.csv(DEGs, file = "DEGs.csv")

# Gene to plot
gene_to_plt <- "KDM6A"
comparison

# Plot normalised counts for a single gene (log-transformed)
plotCounts(dds, gene = gene_to_plt, intgroup = "Reproduction")

# Manually extract normalised counts for gene
gene_to_plt_counts <- normalised_counts[gene_to_plt, ]
# gene_to_plt_counts <- assay(vst)[gene_to_plt, ]

# Create data frame with reproduction information
gene_df <- tibble(
  sample = names(gene_to_plt_counts),
  norm_counts = gene_to_plt_counts,
  group = dds$Reproduction
)

# Convert group to a factor
gene_df$group <- factor(
  x = gene_df$group,
  levels = c("Female_lac_postlac","Female_nonreproductive","Male"),
  labels = c("Females post(lactacting)","Females non-reproductive","Males")
)

# Plot normalised counts coloured by reproduction
set.seed(123)
ggplot(data = gene_df, aes(x = group, y = norm_counts))+
  geom_violin(aes(fill = group), scale = "count", alpha = 0.70, colour = NA)+
  geom_jitter(width = 0.1, fill = "black", shape = 21, size = 3)+
  scale_fill_manual("Group", values = c("#fc8d59","#ffffbf","#91bfdb"))+
  # scale_y_continuous(limits = c(0, NA))+
  ylab("Normalised Transcript Counts")+
  ggtitle(paste("Gene:", gene_to_plt))+
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  )


# ----------------- #
# Temperature model ####
# ----------------- #

# Construct a DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi, coldata_mat, ~ Reproduction + Temperature)

# Estimate range of size factors using Median of ratios normalisation method
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
dds <- estimateSizeFactors(dds)
sort(sizeFactors(dds))
# plot(sizeFactors(dds))

# # Remove outlier samples
# deceased_bat <- c("B680_BL")
# lowSizeFactorRatios <- c("B170710Bba1_BL","B582_BL","B158_BL","B258_BL")
# highSizeFactorRatios <- c("B296_BL","B691_BL","B170713Bba1_BL","B234_BL","B294_BL","B254_BL")
# outlier_samples <- c(deceased_bat, lowSizeFactorRatios, highSizeFactorRatios)
# dds <- dds[, !(colnames(dds) %in% outlier_samples)]
# dds
# 
# # Estimate size factors
# dds <- estimateSizeFactors(dds)
# sort(sizeFactors(dds))

# Print number of samples per group
table(dds[["Temperature"]]); table(dds[["Reproduction"]])

# QC: Gene-level
# Pre-filter: omit genes with zero counts in all samples
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]

# QC: Gene-level
# Pre-filter: omit genes that have less than 10 counts for the smallest group size
# E.g. if minGroupSize=5, five or more samples must have at least 10 counts
minGroupSize <- min(table(dds[["Temperature"]]))
keep <- rowSums(counts(dds) >= 10) >= minGroupSize
dds <- dds[keep, ]
dds

# Total number of counts per sample
# colSums(counts(dds))

# Total number of normalized counts per sample
# colSums(counts(dds, normalized = TRUE))

# Extract normalised counts
normalised_counts <- counts(dds, normalized = TRUE)

# Export normalised counts
# write.csv(normalised_counts, file = "data/normalised_counts.csv", row.names = T)

# Variance Stabilising Transformation (VST)
vst <- varianceStabilizingTransformation(dds, blind = TRUE)
# assay(vst)[1:5, 1:5]

# Export VST counts
# write.csv(assay(vst), file = "data/vst_counts.csv", row.names = T)

# QC: Sample-level
# Principal components analysis
pltPCA(vst, metadata, group = "Temperature", axes = c("PC1","PC2"), size = 5)
pltPCA(vst, metadata, group = "Sex", axes = c("PC1","PC2"), size = 5)
pltPCA(vst, metadata, group = "Reproduction", axes = c("PC1","PC2"), size = 5)
pltPCA(vst, metadata, group = "Country", axes = c("PC1","PC2"), size = 5)

# Run analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Coefficients estimated
resultsNames(dds)

# Levels to compare
comparison <- c("Temperature", "low", "mid") # 20 DEGs
# comparison <- c("Temperature", "low", "high") # 2 DEGs
# comparison <- c("Temperature", "high", "mid") # 3 DEGs

# Extract results table
res <- results(dds, contrast = comparison, alpha = alpha, pAdjustMethod = "BH", test="Wald")
summary(res)

# Shrink log2 fold changes using adaptive shrinkage estimator (ashr) method
res_shrunken <- lfcShrink(dds, contrast = comparison, res = res,type = "ashr")
summary(res_shrunken)

# Plot normalized counts mean versus the log2 foldchanges for all genes tested
# plotMA(res, alpha = alpha, ylim = c(-5,5), colSig = "red", main = "Results")
# plotMA(res_shrunken, alpha = alpha, ylim = c(-5,5), colSig = "red", main = "Results Shrunken")

# Plot log2 fold change versus adjusted p-values
volcano_plt <- as_tibble(res_shrunken) |>
  mutate(.data = _, gene = rownames(res_shrunken)) |>
  mutate(.data = _, DEG = padj < alpha) |>
  drop_na(data = _) |>
  mutate(.data = _, gene_labels = ifelse(DEG == T, gene, "")) |>
  ggplot(data = _)+
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = DEG), size = 3)+
  xlab("log2 fold change")+
  ylab("-log10 adjusted p-value")+
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = gene_labels))
# volcano_plt

# DEGs table
DEGs <- as_tibble(res_shrunken) |>
  # Add column for gene names
  mutate(.data = _, gene = rownames(res_shrunken)) |>
  # Reorder data frame
  dplyr::select(.data = _, gene, everything()) |>
  # Remove columns with NAs
  drop_na(data = _) |>
  # Filter genes with adjusted p-values < alpha
  subset(x = _, padj < alpha) |>
  # Filter genes with absolute values > lfc_threshold
  subset(x = _, abs(log2FoldChange) > lfc_threshold) |>
  # Order by adjusted p-value
  arrange(.data = _, padj)
DEGs
# write.csv(DEGs, file = "DEGs.csv")

# Save DEGs into variables
DEGs1 <- DEGs # Low vs Middle
DEGs2 <- DEGs # Low vs High
DEGs3 <- DEGs # Middle vs High
allDEGs <- rbind(DEGs1, DEGs2, DEGs3)
arrange(allDEGs, log2FoldChange)

# Gene to plot
gene_to_plt <- "IGF1"
comparison

# Plot normalised counts for a single gene (log-transformed)
plotCounts(dds, gene = gene_to_plt, intgroup = "Temperature")

# Manually extract normalised counts for gene
gene_to_plt_counts <- normalised_counts[gene_to_plt, ]
# gene_to_plt_counts <- assay(vst)[gene_to_plt, ]

# Create data frame with metadata information
gene_df <- tibble(
  sample = names(gene_to_plt_counts),
  norm_counts = gene_to_plt_counts,
  group = dds$Temperature,
  reproduction = dds$Reproduction
)

# Convert group to a factor
gene_df$group <- factor(
  x = gene_df$group,
  levels = c("low","mid","high"),
  labels = c("Low","Middle","High")
)

# Plot normalised counts coloured by reproduction
set.seed(123)
ggplot(data = gene_df, aes(x = group, y = norm_counts))+
  geom_violin(aes(fill = group), scale = "count", alpha = 0.70, colour = NA)+
  geom_jitter(width = 0.1, fill = "black", shape = 21, size = 3)+
  scale_fill_manual("Temperatures", values = c("#fc8d59","#ffffbf","#91bfdb"))+
  # scale_y_continuous(limits = c(0, NA))+
  xlab("Temperatures")+
  ylab("Normalised Transcript Counts")+
  ggtitle(paste("Gene:", gene_to_plt))+
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  )

# Plot PCA of VST counts for only DEGs
vst_DEGs_all <- as.data.frame(assay(vst)[allDEGs$gene, ])
# vst_DEGs_males <- dplyr::select(vst_DEGs_all, subset(metadata, Sex == "M")$Sample)
# vst_DEGs_females <- dplyr::select(vst_DEGs_all, subset(metadata, Sex == "F")$Sample)
pca1 <- prcomp(t(vst_DEGs_all))
library(mapmixture)
scatter_plot(as.data.frame(pca1$x), group_ids = metadata$Temperature, type = "labels")
# scatter_plot(as.data.frame(pca1$x), group_ids = subset(metadata, Sex == "M")$Temperature, type = "labels")
# scatter_plot(as.data.frame(pca1$x), group_ids = subset(metadata, Sex == "F")$Temperature, type = "labels")


# ----------------- #
# Gene Enrichment ####
# ----------------- #

# Load libraries
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Explore the AnnotationHub object
ah

# Explore the Data Providers
unique(ah$dataprovider)

# Query AnnotationHub
# https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
species <- "Homo sapiens"
database <- "NCBI"
query_ah <- query(ah, c(species, database))
query_ah

# Extract annotations of interest
# https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
human_ah <- query_ah[["AH111575"]]
human_ah

# Display the columns
columns(human_ah)

# Retrieve gene-level information for all genes tested in model
all_genes <- AnnotationDbi::select(
  x = human_ah,
  keys = rownames(dds),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

# Read in DEGs
DEGs <- read.csv("DEGs.csv") |> distinct()

# Retrieve gene-level information for significant DEGs
sig_genes <- AnnotationDbi::select(
  x = human_ah,
  keys = DEGs$DEG,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

# ENTREZID gene list for background genes and significant genes (DEGs)
all_genes_entrez <- unique(na.exclude(all_genes)[["ENTREZID"]])
sig_genes_entrez <- na.exclude(sig_genes)[["ENTREZID"]]

# Run GO enrichment analysis 
ego <- enrichGO(gene = sig_genes_entrez,
                universe = all_genes_entrez,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db, 
                ont = "ALL",
                pAdjustMethod = "BH", 
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = alpha,
                qvalueCutoff = alpha,
                readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
cluster_summary
write.csv(ego, file = "GO_terms.csv")

# Visualise results
dotplot(ego, showCategory = 50)
emapplot(pairwise_termsim(ego, showCategory = 50))

