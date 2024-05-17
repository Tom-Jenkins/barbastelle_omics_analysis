# ============================ #
#
# Genotype-Environment Association: LFMM2 and RDA
#
# Data source:
# ./data/barb_m50g10maf03_2024.lfmm_imputed.lfmm
# ./sample_coordinates_env.csv
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)

# Load packages
library(LEA)
library(readr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(vcfR)
library(ade4)
library(vegan)
library(stringr)
library(mapmixture)
library(dplyr)

# link to climate data
# https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2Ft
# BIO5 = Max Temperature of Warmest Month
# BIO12 = Annual Precipitation

# Read in genotypes in LFMM format
geno <- data.table::fread("./data/barb_m3_g90_maf03_feb2024.imputed.lfmm")

# Read in sample data and convert env data to a matrix
sample_data <- read_csv("sample_coordinates_env.csv")
env_data <- sample_data[, c("bio5_1981","bio12_1981")]
colnames(env_data) <- c("bio5","bio12")

# Add location factor to sample data
location_ord <- c("Bedfordshire","Dartmoor","Nottinghamshire","Warwickshire","Sussex",
                  "Cazorla","La Rioja","Soria","Teruel","Sabugal")
sample_data$location <- factor(word(sample_data$LOCATION, start = 1, sep = ","), levels = location_ord)

#--------------#
# PCA ####
#--------------#

# Run PCA
pca_geno <- dudi.pca(geno, scale = TRUE, scannf = FALSE, nf = 3)

# Export PCA results as RData file
save(pca_geno, file = "pca_popgen_46230snps.RData")

# Analyse how much percent of genetic variance is explained by each axis
percent <- pca_geno$eig/sum(pca_geno$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)")

# Plot results
# site_ids <- str_remove(sample_data$VCF_ID, "_[A-Z]*[0-9]*")
site_ids <- factor(c(rep("Sussex", 23), "England", rep("Sussex", 12), rep("England", 20), rep("Iberia", 38), "England"), levels = c("England","Sussex","Iberia"))
labels <- str_extract(sample_data$location, "^.{3}")
scatter_plot(
  dataframe = pca_geno$li, group_ids = site_ids, type = "labels", labels = labels,
  colours = c("#91bfdb","#ffffbf","#fc8d59"))+
  theme(
    legend.text = element_text(size = 14),
    legend.position = "top",
    legend.spacing.x = unit(0.5, "cm")
  )

#--------------#
# LFMM2 ####
#--------------#

# Run LEA PCA
pc = pca("./data/barb_m3_g90_maf03_feb2024.imputed.lfmm", scale = TRUE)

# Plot eigenvalues
plot(pc, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")

# Accords with snmf that number of latent factors (K) should be set to 3

# Run LFMM2 analysis
lfmm_mod <- lfmm2(geno, as.matrix(env_data), K = 3, lambda = 1e-5)

# Compute P-values adjusted by genomic inflation factor
lfmm_test_env <- lfmm2.test(lfmm_mod, geno, as.matrix(env_data), linear = T, full = F, genomic.control = T)
lfmm_test_full <- lfmm2.test(lfmm_mod, geno, as.matrix(env_data), linear = T, full = T, genomic.control = T)

# Genomic inflation factor
lfmm_test_env$gif
lfmm_test_full$gif

# Adjusted P-values for full model (one per locus)
lfmm_adj_pval <- p.adjust(lfmm_test_full$pvalues, "BH"); hist(lfmm_adj_pval)

# Adjust P-values for multiple comparisons using Benjamini-Hochberg method
bio5_adj_pval <- p.adjust(lfmm_test_env$pvalues["bio5",], "BH"); hist(bio5_adj_pval)
bio12_adj_pval <- p.adjust(lfmm_test_env$pvalues["bio12",], "BH"); hist(bio12_adj_pval)

#--------------#
# Candidate ddRAD loci
#--------------#

# Expected false discovery rate
alpha <- 0.05

# Data frame of adjusted P-values
df <- data.frame(
  bio5_adj_pval = bio5_adj_pval,
  bio12_adj_pval = bio12_adj_pval,
  full_adj_pval = lfmm_adj_pval
)
head(df)

# Add column indicating if P-value < alpha
df$bio5_sig <- ifelse(df$bio5_adj_pval < alpha, "Signif", "")
df$bio12_sig <- ifelse(df$bio12_adj_pval < alpha, "Signif", "")
df$full_sig <- ifelse(df$full_adj_pval < alpha, "Signif", "")

# Number of candidate loci per variable
table(df$bio5_sig)
table(df$bio12_sig)
table(df$full_sig)

# Row indexes of candidate signif loci
bio5_loci <- which(df$bio5_sig == "Signif")
bio12_loci <- which(df$bio12_sig == "Signif")
full_loci <- which(df$full_sig == "Signif")

# Total number of unique candidate loci
lfmm_loci <- unique(sort.int(c(bio5_loci, bio12_loci, full_loci)))
# lfmm_loci <- full_loci
length(lfmm_loci)

#--------------#
# RDA ####
#--------------#

# Vector of country IDs to represent population structure (K=3)
# env_data$popstructure <- factor(c(rep("Sussex", 23), "England", rep("Sussex", 12), rep("England", 20), rep("Iberia", 38), "England"), levels = c("England","Sussex","Iberia"))
# NOTE: some correlation with BIO5. More outliers detected so currently not included in model.

# Run RDA
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# rda1 <- rda(geno ~ bio5 + bio12 + Condition(popstructure), data = env_data, scale = T)
rda1 <- rda(geno ~ bio5 + bio12, data = env_data, scale = T)
rda1
RsquareAdj(rda1)

# Variance Inflation Factors
vif.cca(rda1)

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "constrained"))

# Screeplot (same number of axes as the number of predictors in the model)
screeplot(rda1)

# Check each canonical axis for significance (**very long run time ~30 mins)
# signif_axis <- anova.cca(rda1, by = "axis", parallel = getOption("mc.cores"))
# signif_axis
# Permutation test for rda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = geno ~ bio5 + bio12, data = env_data, scale = T)
# Df Variance      F Pr(>F)    
# RDA1      1     1910 4.0237  0.001 ***
#   RDA2      1      639 1.3466  0.001 ***
#   Residual 92    43680                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract SNP loadings on constrained axes
rda1_load <- scores(rda1, choices = c(1:3), display = "species")

# Histogram of the loadings on each RDA axis
hist(rda1_load[,1], main="Loadings on RDA1")
hist(rda1_load[,2], main="Loadings on RDA2")

# Function where x is the vector of loadings and z is the number of standard deviations to use
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)    
  x[x < lims[1] | x > lims[2]]
}

# Identify candidate outlier SNPs 
cand1 <- outliers(rda1_load[,1], 3)
cand2 <- outliers(rda1_load[,2], 3)
cand_loci <- unique(names(c(cand1, cand2)))
rda_loci <- as.integer(str_remove(cand_loci, "V"))
length(rda_loci)

#--------------#
# Export Outliers ####
#--------------#

# Outlier loci common among both methods
outlier_loci <- intersect(lfmm_loci, rda_loci)
length(outlier_loci)

# Extract list of candidate loci from vcf
vcf <- read.vcfR("./data/barb_m3_g90_maf03_feb2024.vcf")
vcf@fix[outlier_loci,]
loci_IDs <- paste0("CLocus_", vcf@fix[outlier_loci,][,"CHROM"])
write.table(loci_IDs, file = "outlier_loci.tsv", row.names = F, col.names = F, quote = F)

# Export outlier genotypes in VCF format
write.vcf(vcf[outlier_loci,], "outlier_loci.vcf")

# Export outlier genotypes in lfmm format
geno_outlier <- dplyr::select(geno, all_of(outlier_loci))
write.table(geno_outlier, file = "./data/barb_m3_g90_maf03_feb2024.imputed.lfmm.outliers", row.names = F, col.names = F)

#--------------#
# Outlier RDA ####
#--------------#

# Run RDA with only outlier loci
rda_outlier <- rda(geno_outlier ~ bio5 + bio12, data = env_data, scale = T)

# Total variance explained
RsquareAdj(rda_outlier)

# Create data.frame for plotting
plot_df <- tibble(
  sample = as.factor(sample_data$VCF_ID),
  country = as.factor(site_ids),
  location = sample_data$location,
  bio5 = sample_data$bio5_1981,
  bio12 = sample_data$bio12_1981
)

# Environmnental values per location
plot_df$location

# Vector of colours for bat samples
library(randomcoloR)
library(scales)
sample_cols <- c("#e0f3f8","#74add1","#4575b4","#807dba","#313695",
                 "#ffffbf","#fee090","#e08214","#d73027","#de77ae")
# sample_cols <- distinctColorPalette(n_distinct(plot_df$location), runTsne = FALSE)
scales::show_col(sample_cols)

# Biplot using ordiplot from vegan
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# png("rda_plot.png", width = 10, height = 8, units = "in", res = 600)
ordiplot(rda_outlier, type = "none", choices = c(1,2), scaling = 3) |>
  # points("species", pch = 20, cex = 0.7, col = "grey32") |>
  text("biplot", col = "black", cex = 1) |>
  points("sites", pch = 21, cex = 1.8, col = "black", bg = sample_cols[plot_df$location])

# Add legend to triplot
legend_labels <- str_c(levels(plot_df$location), " ", paste0("(", c(rep("Eng", 5),rep("Spa",4), rep("Por",1)), ")"))
legend("bottomright", legend_labels, bty = "n", col = "black",
       pch = 21, cex = 1, pt.bg = sample_cols)

# Add climate descriptions
mtext_bio5 <- "bio5: maximum temperature of warmest month"
mtext_bio12 <- "bio12: annual precipitation"
mtext(mtext_bio5, adj = 0.01, line = -29, cex = 0.7)
mtext(mtext_bio12, adj = 0.01, line = -30, cex = 0.7)

# Add environmental guides
mtext("Warmer & Drier", adj = 0.05, line = -20, col = "grey70", cex = 1)
mtext("Warmer & Wetter", adj = 0.05, line = -3, col = "grey70", cex = 1)
mtext("Cooler & Wetter", adj = 0.90, line = -3, col = "grey70", cex = 1)

# Add title
title("Redundancy Analysis Biplot")
# dev.off()


#--------------#
# VCF Outlier Genotypes ####
#--------------#

# Extract the outlier genotypes for all individuals
vcf_genotypes <- extract.gt(vcf, return.alleles = TRUE)
vcf_outlier <- vcf_genotypes[vcf@fix[outlier_loci,][,"ID"], ]
nrow(vcf_outlier)
vcf_outlier

# Convert to a data frame
vcf_df <- as.data.frame(t(vcf_outlier))
vcf_df

# Group sample data by highest and lowest BIO5
sample_data %>% select(VCF_ID, bio5_1981) %>% pull(bio5_1981) %>% sort
bio5_grps <- sample_data |>
  select(.data=_, VCF_ID, bio5_1981) |>
  mutate(summer_grp = ifelse(bio5_1981 < 25, "cooler", "warmer")) |>
  pull(summer_grp)
vcf_df$bio5_grps <- bio5_grps

# Group sample data by highest and lowest BIO12
sample_data %>% select(VCF_ID, bio12_1981) %>% pull(bio12_1981) %>% sort
bio12_grps <- sample_data |>
  select(.data=_, VCF_ID, bio12_1981) |>
  mutate(cmi_grp = ifelse(bio12_1981 < 800, "arid", "less_arid")) |>
  pull(cmi_grp)
vcf_df$bio12_grps <- bio12_grps

# Add data for country and sex to vcf data frame
vcf_df$country <- site_ids
vcf_df$sex <- sample_data$SEX
head(vcf_df)

# Extract genotypes for each variable
bio5_geno <- lapply(1:nrow(vcf_outlier), function(x) table(vcf_df[, c(x, which(colnames(vcf_df) == "bio5_grps"))]))
bio12_geno <- lapply(1:nrow(vcf_outlier), function(x) table(vcf_df[, c(x, which(colnames(vcf_df) == "bio12_grps"))]))
country_geno <- lapply(1:nrow(vcf_outlier), function(x) table(vcf_df[, c(x, which(colnames(vcf_df) == "country"))]))
sex_geno <- lapply(1:nrow(vcf_outlier), function(x) table(vcf_df[, c(x, which(colnames(vcf_df) == "sex"))]))

# Add ddRAD loci names to list
names(bio5_geno) <- vcf@fix[outlier_loci,][,"CHROM"]
names(bio12_geno) <- vcf@fix[outlier_loci,][,"CHROM"]
names(country_geno) <- vcf@fix[outlier_loci,][,"CHROM"]
names(sex_geno) <- vcf@fix[outlier_loci,][,"CHROM"]

# Extract table for SNP of interest
snp <- "1267"
bio5_geno[[snp]]; bio12_geno[[snp]]; country_geno[[snp]]; sex_geno[[snp]]

# Extract row index of target outlier loci
target_outlier_loci <- c(MAPKAPK2 = "65203", AUTS2 = "18828", IRX6 = "121099", CUX1 = "116566")
full_outlier_id <- sapply(1:length(target_outlier_loci), function(x) str_subset(colnames(vcf_df), paste0(target_outlier_loci[x], ":[0-9]*")))
target_idx <- match(full_outlier_id, rownames(vcf_genotypes))

# Get row indexes for all outlier loci located in gene regions
ingene_outliers <- c("108280","116566","121099","123012","1267","18828","2620","35311","46103","46844","55654",
                     "55874","62240","63168","64244","65203","90313")
ingene_outliers <- sapply(1:length(ingene_outliers), function(x) str_subset(colnames(vcf_df), paste0(ingene_outliers[x], ":[0-9]*")))
ingene_idx <- match(ingene_outliers, rownames(vcf_genotypes))

#--------------#
# Allele Frequencies ####
#--------------#

# Calculate minor allele frequency for all outliers
(maf <- vcfR::maf(vcf[outlier_loci, ]))

# Group to plot (bio5_grps, bio12_grps, country, sex)
grp <- "bio12_grps"

# Convert to genind object
library(adegenet)
genind_outlier <- df2genind(
  X = vcf_df[ ,1:(ncol(vcf_df)-4)],
  sep = "/",
  ind.names = rownames(vcf_df),
  loc.names = colnames(vcf_df)[1:(length(colnames(vcf_df))-4)],
  pop = vcf_df[[grp]]
)

# Calculate allele frequencies for outlier genotypes per individual
allele_freqs_ind <- makefreq(genind_outlier)

# Calculate average allele frequencies per group
library(poppr)
allele_freq_grp <- as.data.frame(rraf(genind_outlier, by_pop=TRUE, correction = FALSE))
round(allele_freq_grp, digits = 2)
round(allele_freq_grp, digits = 2) |> dplyr::select(contains("1267"))

# Extract country allele frequencies for outlier SNPs in genes
round(dplyr::select(allele_freq_grp, contains(ingene_outliers)), 2)

# Extract country allele frequencies for outlier SNPs in exons
inexon_outliers <- c("35311","46103","63168","121099")
inexon_freq <- dplyr::select(allele_freq_grp, contains(inexon_outliers))
round(inexon_freq, 2)


#--------------#
# Figure 4 ####
#--------------#

# Load libraries
library(ggvenn)
library(tidyr)
library(ggvegan)
library(ggplotify)
library(patchwork)

# List of genes present in both LFMM and RDA
venn_loci <- list(LFMM = lfmm_loci, RDA = rda_loci)

# Plot venn diagram
Fig4A <- ggvenn(
  data = venn_loci,
  fill_color = c("green", "yellow"),
  fill_alpha = 0.1,
  stroke_size = 0.5,
  set_name_size = 3,
  text_size = 3
)
Fig4A

# Data frame for LFMM adj pvals
man_df <- df[,1:3]
head(man_df)

# Add column for smallest adj pval across all three LFMM tests (bio5, bio12 and full)
man_df$adj_pval_min <- apply(man_df, MARGIN = 1, min)
head(man_df)

# Add column value based on row index
man_df$Signif <- ifelse(seq_len(nrow(man_df)) %in% outlier_loci, "SIG", "NONSIG")
head(man_df); nrow(man_df)

# Add column value based on row index
man_df$InGene <- ifelse(seq_len(nrow(man_df)) %in% ingene_idx, "Yes" , "No")
head(man_df); table(man_df$InGene)

# Combine Signif and InGene columns
man_df$Signif_InGene <- case_when(
  man_df$Signif == "NONSIG" ~ "NONSIG",
  man_df$Signif == "SIG" & man_df$InGene == "No" ~ "GEA candidate",
  man_df$Signif == "SIG" & man_df$InGene == "Yes" ~ "GEA candidate (in gene)"
)
head(man_df); table(man_df$Signif_InGene)

# Get position of target gene points
target_pt <- man_df[target_idx,]
target_pt$X <- target_idx
target_pt$Locus <- full_outlier_id
target_pt$Gene <- names(target_outlier_loci)
target_pt

# Manhatten plot
Fig4B <- ggplot(data = man_df)+
  geom_point(aes(x = 1:nrow(man_df), y = -log10(adj_pval_min), colour = Signif_InGene), shape = 16, alpha = 0.7)+
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "yellow")+
  scale_colour_manual(values = c("darkgreen","green2","black"), breaks = c("GEA candidate", "GEA candidate (in gene)"))+
  xlab("SNP")+
  ylab("Adjusted P-value (-log10)")+
  # annotate("label", x = target_pt$X, y = -log10(target_pt$adj_pval_min), label = target_pt$Gene)+
  # theme_minimal()+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = NA),
    legend.spacing.x = unit(0.01, "lines"),
    axis.title = element_text(size = 10),
    panel.grid.major.y = element_line(linewidth = 0.1),
    panel.grid.major.x = element_line(linewidth = 0.1),
    panel.grid.minor.x = element_blank(),
  )
Fig4B

# RDA axes proportion of variance explained
# https://stackoverflow.com/questions/62542609/extracting-proportion-of-variance-explained-from-summaryrda-for-axis-labels
RDA1_percent <- round(summary(rda_outlier)$cont$importance[2,"RDA1"]*100, digits = 1)
RDA2_percent <- round(summary(rda_outlier)$cont$importance[2,"RDA2"]*100, digits = 1)

# Scaling to use
scaling = 3

# Prepare sample points
sites_scores <- fortify(rda_outlier, display = "wa", axes = 1:2, scaling = scaling)
sites_scores$location <- plot_df$location
head(sites_scores)

# Extract data used to plot environmental predictor arrows
arrow_bio5 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][1,]
arrow_bio12 <- ggplot_build(autoplot(rda_outlier, layers = c("sites","biplot"), scaling = scaling))$data[[2]][2,]

# RDA
Fig4C <- ggplot()+
  # Sample points
  geom_point(
    data = sites_scores, aes(x = RDA1, y = RDA2, fill = location),
    shape = 21, colour = "black", size = 4
  )+
  # BIO5 arrow
  annotate("segment", x = 0, xend = arrow_bio5$xend, y = 0, yend = arrow_bio5$yend,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow_bio5$xend-0.35, y = arrow_bio5$yend, label = "bio5", col = "black", size = 5)+
  # BIO12 arrow
  annotate("segment", x = 0, xend = arrow_bio12$xend+0.15, y = 0, yend = arrow_bio12$yend-0.2,
           arrow = arrow(length = unit(0.15, "inches")), colour = "black", linewidth = 0.5)+
  annotate("text", x = arrow_bio12$xend, y = arrow_bio12$yend, label = "bio12", col = "black", size = 5)+
  scale_fill_manual(name = "Location", values = sample_cols)+
  scale_x_continuous(limits = c(-2.5,1.5))+
  scale_y_continuous(limits = c(-3,1.2))+
  xlab(paste0("RDA1 (", RDA1_percent, "%)"))+
  ylab(paste0("RDA2 (", RDA2_percent, "%)"))+
  ggtitle("Redundancy Analysis Biplot")+
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 3)+
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = 3)+
  # Text
  annotate("text", x = -2, y = -2, label = "Warmer & Drier", col = "grey70", size = 4.5)+
  annotate("text", x = -2, y = 1.2, label = "Warmer & Wetter", col = "grey70", size = 4.5)+
  annotate("text", x = 1, y = 1.2, label = "Cooler & Wetter", col = "grey70", size = 4.5)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, colour = "black"),
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12) 
  )
Fig4C


#--------------#
# Figure 4 Export ####
#--------------#

# Layout design
layout <- "
  AACCC
  AACCC
  BBCCC
  BBCCC
"

# Plot layout
plt_list = list(
  Fig4A+ labs(tag = "A"),
  Fig4B+ labs(tag = "B"),
  Fig4C+ labs(tag = "C")
)
Figure4 <- wrap_plots(plt_list, design = layout)
# Figure4

# Export
ggsave(plot = Figure4, filename = "../Figure4.jpg", width = 11, height = 6, dpi = 600)

