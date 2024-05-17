# ============================ #
#
# Barbastelle Differential Gene Expression Analysis
#
#
# Utility R Functions
#
# ============================ #

# ----------------- #
# Principal component analysis
# ----------------- #
pltPCA <- function(matrix, metadata, group, axes = c("PC1","PC2"), colour = NULL,
                   returnData = FALSE, size = 3) {
  
  # Perform PCA
  pca <- prcomp(t(SummarizedExperiment::assay(matrix)))
  
  # Subset metadata data frame
  samples <- colnames(matrix)
  metadata <- subset(metadata, subset = metadata$Sample %in% colnames(vst))
  
  # Data frame to merge sample metadata and PCs
  df <- cbind(metadata, pca$x)
  
  # Scatter plot using ggplot2
  plt <- ggplot(data = df)+
    geom_point(
      aes(x = !!as.name(axes[1]), y = !!as.name(axes[2]), color = factor(!!as.name(group))),
      size = size,
    )+
    guides(colour = guide_legend(title = group))
  
  # Return plot
  if (!returnData) {
    return(plt)
  }  
  
  # Return data
  if (returnData) {
    return_df <- df[, axes]
    return(return_df)
  }
}

# ----------------- #
# Hierarchical clustering
# ----------------- #
pheatmap2 <- function(matrix, coldata, group, legend_title = "", ...) {
  
  # Correlation matrix
  cor_mat <- cor(matrix)
  
  # Title of legend
  if (legend_title == "") {
    legend_title = group
  }
  
  # Annotation data frame
  annotate_df <- as.data.frame(coldata[, group])
  colnames(annotate_df) <- legend_title
  
  # pheatmap function
  pheatmap(
    mat = cor_mat,
    annotation_col = annotate_df,
    ...
  )
}



