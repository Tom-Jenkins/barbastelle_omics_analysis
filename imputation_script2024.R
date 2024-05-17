# ============================ #
#
# Genotype-Environment Association: Imputation
#
# Data source:
# ./data/barb_m50g10maf03_2024.vcf
#
# ============================ #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(LEA)
library(data.table)

# Convert vcf to lfmm format
# vcf2lfmm("./data/barb_m50g10maf03_2024.vcf", "./data/barb_m50g10maf03_2024.lfmm")
vcf2lfmm("./data/barb_m3_g90_maf03_feb2024.vcf", "./data/barb_m3_g90_maf03_feb2024.vcf")

#--------------#
# SNMF
#--------------#

# K = 3 determined from the snmf analysis by Razgour et al. 2023
# https://doi.org/10.1111/1365-2664.14540

# Perform snmf analysis
project <- snmf(
  input.file = "./data/barb_m3_g90_maf03_feb2024.vcf",
  K = 3,
  project = "new",
  repetitions = 10,
  entropy = TRUE,
  ploidy = 2
)

#--------------#
# Imputation
#--------------#

# Select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project, K = 3))

# Impute the missing genotypes
impute(
  object = project,
  input.file = "./data/barb_m3_g90_maf03_feb2024.lfmm",
  method = "mode",
  K = 3,
  run = best
)

#--------------#
# Remove monomorphic loci
# Note: this step was required for downstream analyses to work
#--------------#

# Read in imputed file
imputed_loci <- read.lfmm("./data/barb_m3_g90_maf03_feb2024.lfmm_imputed.lfmm")

# Get rows where all values are identical (monomorphic SNPs)
row_uniq_vals <- apply(imputed_loci, 2, data.table::uniqueN)
row_polymorphic <- which(row_uniq_vals > 1)

# Subset matrix for polymorphic loci
imputed_loci_poly <- imputed_loci[ ,row_polymorphic]

#--------------#
# Export imputed and polymorphic loci
#--------------#

# Export matrix in lfmm format
write.lfmm(imputed_loci_poly, "./data/barb_m3_g90_maf03_feb2024.imputed.lfmm")

# Export in geno format
lfmm2geno("./data/barb_m3_g90_maf03_feb2024.imputed.lfmm")
