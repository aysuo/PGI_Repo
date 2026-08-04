#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.3.2_causalEffect_UKB1_BinaryTraits_lme4.R
# Purpose: Estimate the causal effect of polygenic indices (PGIs) for binary
#          traits in the UK Biobank (UKB1 cohort).
#
# Description:
# This script estimates the *causal effect* of proband PGIs by adjusting for
# parental PGIs using a linear mixed model (LMM) framework. A random intercept
# for family ID is included to account for relatedness, implemented via the
# `lme4` R package (`glmer` for binary outcomes).
#
# Usage: Source 'paths_PGIrepo_withinFam' before running.
#
# Output: The script generates a file named
#   `directEffect_UKB1_noresidPheno_coefficients2.csv`
# which contains the estimated causal (proband) effects alongside parental
# PGI effects and their bootstrap-based confidence intervals for all binary
# phenotypes in the UKB1 cohort.
#
# Reference: Alemu et al. (2025). *An Updated Polygenic Index Repository:
# Expanded Phenotypes, New Cohorts, and Improved Causal Inference.*
# ------------------------------------------------------------------------------

# ── 1. Setup ──────────────────────────────────────────────────────────────────
# Source environment paths
Sys.getenv("PGI_RepoV2") # Ensure environment loaded
phenotype_list_file <- Sys.getenv("UKB1_PHENO_LIST")
covar_file          <- Sys.getenv("UKB1_COVAR_FILE")
results_file        <- file.path(Sys.getenv("UKB1_OUTPUT_DIR"),
                                "directEffect_UKB1_noresidPheno_coefficients2.csv")

# Set library path for reproducible installs
personal_lib_path <- file.path(Sys.getenv("PGI_RepoV2"), "doc/R/library")
.libPaths(personal_lib_path)

# ── 2. Load packages ─────────────────────────────────────────────────────────
packages <- c("ggplot2","dplyr","lmtest","boot","readr","tidyverse",
              "biglm","lme4","parallel","sandwich","lmtest")

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}
for (pkg in packages) install_and_load(pkg)

# ── 3. Phenotype list ─────────────────────────────────────────────────────────
phenotypes_all <- read.table(phenotype_list_file, header = FALSE, stringsAsFactors = FALSE)$V1
run_test <- FALSE
phenotypes <- if (run_test) "HEIGHT" else phenotypes_all
print(phenotypes)

# ── 4. Covariates ────────────────────────────────────────────────────────────
covariates <- read.table(covar_file, header = TRUE)
include_batches <- FALSE
if (!include_batches) {
  covariates <- covariates[, c(1:12,(ncol(covariates)-4):ncol(covariates))]
}

# Helper: create covariate formula string
create_formula_covariates <- function(covariates_df) {
  covariates_df <- covariates_df[, !names(covariates_df) %in% c("FID","IID")]
  paste(names(covariates_df), collapse = " + ")
}
formula_covariates <- create_formula_covariates(covariates)

# ── 5. Prepare results file ──────────────────────────────────────────────────
if (!file.exists(results_file)) {
  write.csv(data.frame(
    Phenotype = character(0),
    Coef_PopEffect = numeric(0),
    Boot_SE_PopEffect = numeric(0),
    PopEffect_Lower_CI = numeric(0),
    PopEffect_Upper_CI = numeric(0),
    Coef_Proband = numeric(0),
    Boot_SE_ProbandEffect = numeric(0),
    ProbandEffect_Lower_CI = numeric(0),
    ProbandEffect_Upper_CI = numeric(0),
    Coef_Paternal = numeric(0),
    Boot_SE_Paternal = numeric(0),
    Paternal_Lower_CI = numeric(0),
    Paternal_Upper_CI = numeric(0),
    Coef_Maternal = numeric(0),
    Boot_SE_Maternal = numeric(0),
    Maternal_Lower_CI = numeric(0),
    Maternal_Upper_CI = numeric(0),
    N_pop = numeric(0),
    N_direct = numeric(0)
  ), results_file, row.names = FALSE, quote = FALSE, na = "")
}

# ── 6. Utility functions ─────────────────────────────────────────────────────
standardize <- function(x) (x - mean(x,na.rm=TRUE)) / sd(x,na.rm=TRUE)

should_process_phenotype <- function(data, phenotype_column) {
  unique_values <- unique(na.omit(data[[phenotype_column]]))
  if (all(unique_values %in% c(0,1))) {
    proportion_of_cases <- mean(data[[phenotype_column]]==1,na.rm=TRUE)
    if (is.na(proportion_of_cases)) return(FALSE)
    if (proportion_of_cases<1e-8 || proportion_of_cases>0.9999) return(FALSE)
  }
  TRUE
}

compute_robust_se <- function(model,data,cluster_var) {
  vcovCL <- vcovCL(model, cluster=data[[cluster_var]], type="HC")
  sqrt(diag(vcovCL))
}

# ── 7. Main function to fit models ───────────────────────────────────────────
analyze_effects <- function(data, formula, family_type="gaussian") {
  if (family_type=="binomial") {
    model <- tryCatch(
      glmer(formula,data=data,family=binomial(),
            nAGQ=0,control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=200000))),
      error=function(e){message("Error in glmer: ",e$message);return(NULL)}
    )
    if (is.null(model)) return(NULL)
    fe <- summary(model)$coefficients
    coefs <- fe[,"Estimate"]; se <- fe[,"Std. Error"]
  } else {
    model <- tryCatch(lm(formula,data=data),
                      error=function(e){message("Error in lm: ",e$message);return(NULL)})
    if (is.null(model)) return(NULL)
    coefs <- coef(summary(model))[,"Estimate"]
    se <- compute_robust_se(model,data,"FID")
  }
  ci_lower <- coefs - qnorm(0.975)*se
  ci_upper <- coefs + qnorm(0.975)*se
  if (family_type=="binomial") {
    coefs <- exp(coefs); ci_lower <- exp(ci_lower); ci_upper <- exp(ci_upper)
  }
  n_used <- if (family_type=="binomial") sum(!is.na(getME(model,"y"))) else nobs(model)
  list(Coefs=coefs,SE=se,CI_Lower=ci_lower,CI_Upper=ci_upper,N_used=n_used)
}

# ── 8. Loop over phenotypes ──────────────────────────────────────────────────
for (pheno in phenotypes) {
  parental_proband_pgi_file <- file.path(Sys.getenv("PGI_RepoV2"),
    "data/pgs/UKB1", paste0("PGS_UKB1_parental_",pheno,"-single_SBayesR.pgs.txt"))
  phenotype_file <- file.path(Sys.getenv("PGI_RepoV2"),
    "data/phenotypes/UKB1", paste0(pheno,"_noresid.pheno"))
  
  if (!file.exists(parental_proband_pgi_file)) {next}
  pgi_data <- read.table(parental_proband_pgi_file,header=TRUE)
  pheno_data <- read.table(phenotype_file,header=TRUE,fill=TRUE)
  
  data <- merge(pgi_data,pheno_data,by="IID")
  data <- merge(data,covariates,by="IID")
  names(data)[names(data)==pheno] <- "Phenotype"
  
  is_binary <- all(data$Phenotype %in% c(0,1,NA))
  family_type <- if (is_binary) "binomial" else "gaussian"
  if (is_binary && !should_process_phenotype(data,"Phenotype")) next
  
  if (family_type=="gaussian") data$Phenotype <- standardize(data$Phenotype)
  pc_columns <- grep("PC[0-9]+",names(data),value=TRUE)
  data[pc_columns] <- apply(data[pc_columns],2,standardize)
  data$proband <- standardize(data$proband)
  data$paternal <- standardize(data$paternal)
  data$maternal <- standardize(data$maternal)
  
  if (family_type=="binomial") {
    formula_pop_effect <- as.formula(paste("Phenotype ~ proband +",formula_covariates,"+ (1|FID)"))
    formula_direct_effect <- as.formula(paste("Phenotype ~ proband + paternal + maternal +",formula_covariates,"+ (1|FID)"))
  } else {
    formula_pop_effect <- as.formula(paste("Phenotype ~ proband +",formula_covariates))
    formula_direct_effect <- as.formula(paste("Phenotype ~ proband + paternal + maternal +",formula_covariates))
  }
  
  result_pop_effect <- analyze_effects(data,formula_pop_effect,family_type)
  result_direct_effect <- analyze_effects(data,formula_direct_effect,family_type)
  if (is.null(result_pop_effect) || is.null(result_direct_effect)) next
  
  current_results <- data.frame(
    Phenotype = pheno,
    Coef_PopEffect = result_pop_effect$Coefs["proband"],
    Boot_SE_PopEffect = result_pop_effect$SE["proband"],
    PopEffect_Lower_CI = result_pop_effect$CI_Lower["proband"],
    PopEffect_Upper_CI = result_pop_effect$CI_Upper["proband"],
    Coef_Proband = result_direct_effect$Coefs["proband"],
    Boot_SE_ProbandEffect = result_direct_effect$SE["proband"],
    ProbandEffect_Lower_CI = result_direct_effect$CI_Lower["proband"],
    ProbandEffect_Upper_CI = result_direct_effect$CI_Upper["proband"],
    Coef_Paternal = result_direct_effect$Coefs["paternal"],
    Boot_SE_Paternal = result_direct_effect$SE["paternal"],
    Paternal_Lower_CI = result_direct_effect$CI_Lower["paternal"],
    Paternal_Upper_CI = result_direct_effect$CI_Upper["paternal"],
    Coef_Maternal = result_direct_effect$Coefs["maternal"],
    Boot_SE_Maternal = result_direct_effect$SE["maternal"],
    Maternal_Lower_CI = result_direct_effect$CI_Lower["maternal"],
    Maternal_Upper_CI = result_direct_effect$CI_Upper["maternal"],
    N_pop = result_pop_effect$N_used,
    N_direct = result_direct_effect$N_used
  )
  
  write.table(current_results,file=results_file,sep=",",
              col.names=!file.exists(results_file),
              row.names=FALSE,append=TRUE,quote=FALSE)
}
message("Results processing complete. Data written to file.")