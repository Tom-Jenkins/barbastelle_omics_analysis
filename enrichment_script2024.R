# ============================ #
#
# GEA Outlier Gene Enrichment Analysis
#
# Data source:
# GEA_genes.txt
# Myotis_genes.txt
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load libraries
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

# Read in gene lists
gea_genes <- read.table("GEA_table.tsv", header = F) |> dplyr::select(last_col()) |> setnames("Gene")
myotis_genes <- read.delim("Myotis_genes.txt", header = F, col.names = c("Gene"))

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
myotis_genes <- AnnotationDbi::select(
  x = human_ah,
  keys = myotis_genes$Gene,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

# Retrieve gene-level information for significant DEGs
gea_genes <- AnnotationDbi::select(
  x = human_ah,
  keys = unique(gea_genes$Gene),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

# ENTREZID gene list for background genes and GEA genes
all_genes_entrez <- unique(na.exclude(myotis_genes)[["ENTREZID"]])
gea_genes_entrez <- na.exclude(gea_genes)[["ENTREZID"]]

# Run GO enrichment analysis 
ego <- enrichGO(gene = gea_genes_entrez,
                universe = all_genes_entrez,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db, 
                ont = "ALL",
                pAdjustMethod = "BH", 
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
cluster_summary
# write.csv(ego, file = "GO_terms.csv")

# Visualise results
dotplot(ego, showCategory = 50)
emapplot(pairwise_termsim(ego, showCategory = 50))

