#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.3.3_populationAssociation_UKB3_QuantitativeTraits_lm.R
# Purpose: Estimate the population association of polygenic indices (PGIs) for
#          quantitative traits in the UK Biobank (UKB3 cohort).
#
# Description:
# This script estimates the *population association* of proband PGIs, defined as
# the coefficient on the individual's PGI in a regression of the phenotype on that
# PGI without controls for parental PGIs.
#
# Analyses are conducted in the UKB3 sample — the third partition of the UK Biobank
# consisting of unrelated individuals. Models are estimated using ordinary least
# squares (OLS) regression and include 20 genetic principal components (PCs) as
# covariates to adjust for population structure.
#
# Usage: Source 'paths_PGIrepo_withinFam' before running.
#
# Output: The script generates a file named
#   `popEffect_coefficients.txt`
# which contains the estimated standardized PGI coefficient (population association),
# its standard error, and the sample size for each quantitative phenotype in UKB3.
#
# Reference: Alemu et al. (2025). *An Updated Polygenic Index Repository:
# Expanded Phenotypes, New Cohorts, and Improved Causal Inference.*
# ------------------------------------------------------------------------------

# Source external path definitions
source("paths_PGIrepo_withinFam")

# Set personal R library path
.libPaths(personal_lib_path)

# List of required packages
packages <- c("dplyr", "readr", "lmtest", "boot", "tidyverse", "biglm", "sandwich")

# Install and load packages if not already available
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
  }
  library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
}

# Read the phenotype list
phenotypes <- read.table(ukb3_pheno_list_file, header = FALSE, stringsAsFactors = FALSE)$V1

# Read and process the PCs file
pcs <- read.table(ukb3_covariate_file, header = TRUE)[, 1:22]
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:20))

# Convert PC columns to numeric safely
for (i in 3:22) {
  pcs[, i] <- as.numeric(as.character(pcs[, i]))
  if (any(is.na(pcs[, i]))) {
    warning(paste("NAs introduced in PC column:", colnames(pcs)[i]))
  }
}

# Initialize output data frame
results_df <- data.frame(Phenotype = character(),
                         Coef_Proband = numeric(),
                         SE_Proband = numeric(),
                         N_Pop = integer(),
                         stringsAsFactors = FALSE)

# Function to clean phenotype name (remove trailing numbers)
remove_trailing_numbers <- function(name) {
  sub("[0-9]+$", "", name)
}

# Define model formula for population association
formula_pop_effect <- as.formula("phenotype ~ proband + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                                  PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")

# Loop through phenotypes
for (pheno in phenotypes) {

  pheno_modified <- remove_trailing_numbers(pheno)

  # Define file paths
  pgi_file <- file.path(ukb3_pgi_dir, paste0("PGS_UKB3_", pheno_modified, "-single_SBayesR.txt"))
  phenotype_file <- file.path(ukb3_pheno_dir, paste0(pheno_modified, ".pheno"))

  # Check file existence
  if (!file.exists(phenotype_file)) {
    cat("Phenotype file missing for", pheno, "- skipping.\n")
    next
  }
  if (!file.exists(pgi_file)) {
    cat("PGI file missing for", pheno, "- skipping.\n")
    next
  }

  # Load data
  pgi_data <- readr::read_delim(pgi_file, delim = "\t", col_types = readr::cols())
  phenotype_data <- read.table(phenotype_file, header = TRUE, sep = " ")
  phenotype_data$phenotype <- as.numeric(phenotype_data$phenotype)

  # Rename last column to "proband" (PGI column)
  names(pgi_data)[ncol(pgi_data)] <- "proband"

  # Standardize PGI and phenotype
  pgi_data$proband <- scale(pgi_data$proband)
  phenotype_data$phenotype <- scale(phenotype_data$phenotype)

  # Merge data
  merged_data <- merge(pgi_data, phenotype_data, by = c("FID", "IID"))
  merged_data <- merge(merged_data, pcs, by = c("FID", "IID"))

  # Fit linear model
  fit <- lm(formula_pop_effect, data = merged_data)

  # Extract estimates
  coef_proband <- coef(fit)["proband"]
  se_proband <- summary(fit)$coefficients["proband", "Std. Error"]
  n_obs <- nobs(fit)

  # Store results
  results_df <- rbind(results_df,
                      data.frame(Phenotype = pheno_modified,
                                 Coef_Proband = coef_proband,
                                 SE_Proband = se_proband,
                                 N_Pop = n_obs))
  
  cat("Processed:", pheno, "\n")
}

# Save results
write.table(results_df, file = file.path(ukb3_output_dir, "popEffect_coefficients.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Results saved to: ", file.path(ukb3_output_dir, "popEffect_coefficients.txt"), "\n")