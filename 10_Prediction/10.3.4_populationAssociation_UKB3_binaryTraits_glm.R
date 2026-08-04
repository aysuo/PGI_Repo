#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.3.4_populationAssociation_UKB3_BinaryTraits_glm.R
# Purpose: Estimate the population association of polygenic indices (PGIs) for
#          binary traits in the UK Biobank (UKB3 cohort).
#
# Description:
# This script estimates the *population association* of proband PGIs, defined as
# the coefficient on the individual's PGI in a regression of the phenotype on that
# PGI without controls for parental PGIs. For binary traits, models are estimated
# via logistic regression (glm with binomial link).
#
# Analyses are conducted in the UKB3 sample — the third partition of the UK Biobank
# consisting of unrelated individuals.
#
# Usage: Source 'paths_PGIrepo_withinFam' before running.
#
# Output: The script generates a file named
#   `popEffect_noresidPheno_coefficients3.txt`
# which contains the estimated odds ratios for PGIs, their standard errors, and
# 95% confidence intervals for each binary phenotype in UKB3.
#
# Reference: Alemu et al. (2025). *An Updated Polygenic Index Repository:
# Expanded Phenotypes, New Cohorts, and Improved Causal Inference.*
# ------------------------------------------------------------------------------

# Load custom library path
personal_lib_path <- file.path(Sys.getenv("PGI_RepoV2"), "doc/R/library")
.libPaths(personal_lib_path)

# Install and load necessary packages
packages <- c("dplyr", "readr", "lmtest", "boot", "tidyverse", "biglm", "sandwich", "stringr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
  }
  library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
}

# Read list of phenotypes
phenotypes_all <- read.table(file.path(Sys.getenv("PGI_RepoV2"), "doc/PhenoList/version_single_UKB3"),
                             header = FALSE, stringsAsFactors = FALSE)$V1

# Use all phenotypes unless testing
testing <- FALSE
phenotypes <- if (testing) "XX" else phenotypes_all
print(phenotypes)

# Load covariates
covar_file <- file.path(Sys.getenv("UKB3_covar"), "UKB3.covar")
covariates <- read.table(covar_file, header = TRUE)

# Optionally remove batch covariates
include_batches <- FALSE
if (!include_batches) {
  covariates <- covariates[, c(1:22, (ncol(covariates)-4):ncol(covariates))]
}

# Build regression formula
create_formula <- function() {
  covariate_names <- colnames(covariates)[-c(1, 2)]
  covariates_formula <- paste(covariate_names, collapse = " + ")
  as.formula(paste("phenotype ~ proband +", covariates_formula))
}

remove_trailing_numbers <- function(name) {
  sub("[0-9]+$", "", name)
}

# Output file
results_file <- file.path(Sys.getenv("PGI_RepoV2"),
  "processed/sumStats/within_family_related/UKB3/popEffect_noresidPheno_coefficients3.txt")
if (file.exists(results_file)) file.remove(results_file)

# Loop over phenotypes
for (pheno in phenotypes) {
  pheno_modified <- remove_trailing_numbers(pheno)

  # Set PGI path
  if (pheno_modified == "ALZnoproxy") {
    pgi_file <- file.path(Sys.getenv("PGI_Path_UKB3"), "PGS_UKB3_ALZ-single_SBayesR.txt")
  } else if (pheno_modified == "CAD") {
    pgi_file <- file.path(Sys.getenv("PGI_Path_UKB3"), "PGS_UKB3_CAD-single_SBayesR.txt")
  } else {
    pgi_file <- file.path(Sys.getenv("PGI_Path_UKB3"), paste0("PGS_UKB3_", pheno_modified, "-single_SBayesR.txt"))
  }

  # Set phenotype file path
  phenotype_file <- file.path(Sys.getenv("Phenotype_Path_UKB3"),
                               paste0(pheno_modified, "_noresid.pheno"))

  # Skip if missing
  if (!file.exists(phenotype_file)) {
    cat("Phenotype file for", pheno, "not found. Skipping.\n")
    next
  }
  if (!file.exists(pgi_file)) {
    cat("PGI file for", pheno, "not found. Skipping.\n")
    next
  }

  # Load data
  pgi_data <- read_delim(pgi_file, delim = "\t", col_types = cols())
  phenotype_data <- read.table(phenotype_file, header = TRUE)
  if (pheno == "CAD" && "SOFTCAD" %in% colnames(phenotype_data)) {
    colnames(phenotype_data)[colnames(phenotype_data) == "SOFTCAD"] <- "CAD"
  }
  phenotype_column <- colnames(phenotype_data)[3]
  phenotype_data$phenotype <- as.numeric(phenotype_data[[phenotype_column]])

  # Identify binary traits
  non_na_vals <- na.omit(phenotype_data$phenotype)
  is_binary <- all(non_na_vals %in% c(0, 1))
  if (!is_binary) next

  # Standardize PGI
  names(pgi_data)[ncol(pgi_data)] <- "proband"
  pgi_data$proband <- scale(pgi_data$proband)

  # Merge datasets
  merged_data <- merge(merge(pgi_data, phenotype_data, by = c("FID", "IID")), covariates, by = c("FID", "IID"))

  # Fit logistic regression
  formula_pop_effect <- create_formula()
  fit <- glm(formula_pop_effect, data = merged_data, family = binomial(), control = glm.control(maxit = 50))
  n_used <- nobs(fit)

  # Extract effect
  coef_proband <- exp(coef(fit)["proband"])
  se_proband <- summary(fit)$coefficients["proband", "Std. Error"]
  ci_log_odds <- confint(fit)["proband", ]
  lower_ci <- exp(ci_log_odds[1])
  upper_ci <- exp(ci_log_odds[2])

  # Save to file
  results_row <- data.frame(
    Phenotype = pheno,
    Coef_Proband = coef_proband,
    SE_Proband = se_proband,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    N_Pop = n_used
  )
  write.table(results_row, file = results_file, append = TRUE, sep = "\t",
              row.names = FALSE, quote = FALSE, col.names = !file.exists(results_file))

  cat("Processed phenotype:", pheno, "\n")
}

cat("All results saved to", results_file, "\n")