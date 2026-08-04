#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.4.1_causalEffect_WLS_BinaryTraits_lme4.R
# Purpose: Estimate the causal effect of polygenic indices (PGIs) for binary
#          traits in the Wisconsin Longitudinal Study (WLS).
#
# Description:
#   This script estimates the *causal effect* of proband PGIs by adjusting for
#   parental PGIs in a generalized linear mixed model (GLMM) framework with a
#   family-level random intercept (implemented via `lme4::glmer` with a
#   binomial link).
#
#   For each binary phenotype:
#     - Population association model:
#         phenotype ~ proband + PCs + age/sex terms + (1 | familyID)
#     - Within-family (causal) model:
#         phenotype ~ proband + parental + PCs + age/sex terms + (1 | familyID)
#
#   The script returns log-odds coefficients, their standard errors, and
#   exponentiated odds ratios for proband and parental PGIs, plus sample size.
#
# Usage:
#   1) Ensure the following environment variables are set (e.g. via
#      `paths_PGIrepo_withinFam` bash file):
#         PGI_RepoV2        : base directory of the project
#         WLS_PHENO_LIST    : list of WLS phenotypes (one per line)
#         WLS_AGESEX_FILE   : agesex covariate file
#         WLS_PC_FILE       : principal components file
#         WLS_PGI_DIR       : directory with parental+proband PGI files
#         WLS_PHENO_NORESD_DIR : directory with binary (noresid) phenotypes
#         WLS_OUTPUT_DIR    : directory for WLS causal-effect outputs
#
#   2) Run:
#        Rscript 10.4.1_causalEffect_WLS_BinaryTraits_lme4.R
#
# Output:
#   - <WLS_OUTPUT_DIR>/directEffect_WLS_noresidPheno_coefficients.csv
#     Columns: Phenotype, Coef_PopEffect, SE_PopEffect, OddsRatio_PopEffect,
#              Coef_Proband, SE_Proband, OddsRatio_Proband,
#              Coef_Parental, SE_Parental, OddsRatio_Parental, N
#
# Reference:
#   Alemu et al. (2025). "An Updated Polygenic Index Repository:
#   Expanded Phenotypes, New Cohorts, and Improved Causal Inference."
# ------------------------------------------------------------------------------

# ── 1. Setup ──────────────────────────────────────────────────────────────────

# Ensure project root is defined
PGI_RepoV2 <- Sys.getenv("PGI_RepoV2")
if (PGI_RepoV2 == "") {
  stop("PGI_RepoV2 is not defined. Please set it in your environment before running.")
}

# Input paths from environment (set in paths_PGIrepo_withinFam bash file)
phenotype_list_file <- Sys.getenv("WLS_PHENO_LIST")
agesex_file         <- Sys.getenv("WLS_AGESEX_FILE")
pc_file             <- Sys.getenv("WLS_PC_FILE")
pgi_dir             <- Sys.getenv("WLS_PGI_DIR")
pheno_n_resid_dir   <- Sys.getenv("WLS_PHENO_NORESD_DIR")
wls_output_dir      <- Sys.getenv("WLS_OUTPUT_DIR")

if (any(c(phenotype_list_file, agesex_file, pc_file, pgi_dir,
          pheno_n_resid_dir, wls_output_dir) == "")) {
  stop("One or more required WLS path environment variables are not set.\n",
       "Required: WLS_PHENO_LIST, WLS_AGESEX_FILE, WLS_PC_FILE, ",
       "WLS_PGI_DIR, WLS_PHENO_NORESD_DIR, WLS_OUTPUT_DIR.")
}

results_file <- file.path(
  wls_output_dir,
  "directEffect_WLS_noresidPheno_coefficients.csv"
)

# Set library path for reproducible installs
personal_lib_path <- file.path(PGI_RepoV2, "doc/R/library")
.libPaths(c(personal_lib_path, .libPaths()))

# ── 2. Load packages ─────────────────────────────────────────────────────────

packages <- c(
  "ggplot2", "dplyr", "lmtest", "boot", "readr", "tidyverse",
  "broom", "sandwich", "lme4"
)

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}
for (pkg in packages) install_and_load(pkg)

# ── 3. Read phenotype list and covariates ────────────────────────────────────

# Phenotype names to analyze
phenotypes <- readr::read_lines(phenotype_list_file)
phenotypes <- phenotypes[phenotypes != ""]

message("Number of phenotypes in WLS list: ", length(phenotypes))

# Agesex covariates
agesex_data <- read.table(
  agesex_file,
  header = TRUE,
  fill   = TRUE,
  na.strings = ""
)

# IID clean-up (drop suffix after underscore)
agesex_data$IID <- sapply(strsplit(as.character(agesex_data$IID), "_"), `[`, 1)

# Remove duplicates (keep first occurrence)
agesex_data <- agesex_data %>%
  dplyr::group_by(IID) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# Convert sex: M/F → 0/1 (M = 0, F = 1)
agesex_data$sex <- as.numeric(agesex_data$sex == "F")

# Derived covariates
agesex_data$age2    <- agesex_data$age^2
agesex_data$sexage  <- agesex_data$sex * agesex_data$age
agesex_data$sexage2 <- agesex_data$sex * agesex_data$age2

message("Agesex covariates loaded. Example rows:")
print(utils::head(agesex_data))

# PCs
pcs <- read.table(pc_file, header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))
pcs$true_id <- sapply(strsplit(as.character(pcs$IID), "_"), `[`, 1)
pcs <- pcs %>%
  dplyr::distinct(true_id, .keep_all = TRUE)

message("PCs loaded. Example rows:")
print(utils::head(pcs))

# ── 4. Define phenotype categories (for reference / consistency) ─────────────

phenotype_categories <- list(
  "Anthropometric" = c("BMI", "HEIGHT"),
  "Blood Biomarkers" = c("BL_HDL", "BL_LDL", "BL_nonHDL", "BL_CHOL", "BL_TRYG",
                         "BPdia", "BPpulse", "BPsys"),
  "Cognition & Education" = c("CP", "EA", "ALZnoproxy", "ALZ"),
  "Personality & Wellbeing" = c("ADVENTURE", "EXTRA", "FAMSAT", "FRIENDSAT",
                                "MORNING", "NEURO", "OPEN", "ACTIVITY",
                                "RELIGATT", "RISK", "SWB"),
  "Health" = c("ALLERGYPOLLEN", "ASTHMA", "ASTECZRHI", "BRCA", "COPD",
               "HARDCAD", "HAYFEVER", "IBD", "MIGRAINE", "NEARSIGHTED",
               "PRCA", "SELFHEALTH", "T2D"),
  "Fertility & Sexual Development" = c("AFB", "MENARCHE", "AFS", "NEBwomen"),
  "Psychiatric Conditions" = c("ADHD", "ANOREX", "ASD", "BIPOLAR", "DEP",
                               "INSOMNIA", "SCZ"),
  "Substance Use" = c("AUDIT", "ASI", "CANNABIS", "CPD", "DPW", "EVERSMOKE", "SMCESS")
)

# ── 5. Initialize results container ──────────────────────────────────────────

results_list <- list()

# ── 6. Main loop over phenotypes (binary traits only) ────────────────────────

for (pheno in phenotypes) {
  message("────────────────────────────────────────────")
  message("Processing phenotype: ", pheno)

  # 6.1 Phenotype file (binary: noresid directory)
  phenotype_file <- file.path(pheno_n_resid_dir, paste0(pheno, ".pheno"))
  if (!file.exists(phenotype_file)) {
    message("  ✗ Binary (noresid) phenotype file not found for ", pheno,
            ". Skipping (binary script).")
    next
  }

  # 6.2 PGI file (proband + parental)
  if (pheno == "HARDCAD") {
    pgi_file <- file.path(
      pgi_dir,
      "PGS_WLS_parental_CAD-single_SBayesR.pgs.txt"
    )
  } else {
    pgi_file <- file.path(
      pgi_dir,
      paste0("PGS_WLS_parental_", pheno, "-single_SBayesR.pgs.txt")
    )
  }

  if (!file.exists(pgi_file)) {
    message("  ✗ PGI file not found for ", pheno, ": ", pgi_file, ". Skipping.")
    next
  }

  # 6.3 Load PGI data and scale
  parental_proband_pgi <- read.table(pgi_file, header = TRUE)

  parental_proband_pgi$true_id <- sapply(
    strsplit(as.character(parental_proband_pgi$IID), "_"), `[`, 1
  )

  parental_proband_pgi <- parental_proband_pgi %>%
    dplyr::mutate(
      proband  = as.numeric(scale(proband)),
      parental = as.numeric(scale(parental))
    )

  if (any(duplicated(parental_proband_pgi$true_id))) {
    message("  ⚠ Duplicates in true_id of PGI data for ", pheno,
            " (keeping all rows; family clustering handled in model).")
  }

  # 6.4 Load phenotype data
  pheno_dat <- read.table(phenotype_file, header = TRUE, sep = ",")
  if ("pheno" %in% names(pheno_dat)) {
    names(pheno_dat)[names(pheno_dat) == "pheno"] <- "phenotype"
  }

  pheno_dat$IID <- as.character(pheno_dat$id)
  pheno_dat <- pheno_dat %>%
    dplyr::distinct(IID, .keep_all = TRUE)

  if (!all(pheno_dat$phenotype %in% c(0, 1, NA))) {
    message("  ⚠ Phenotype values not strictly 0/1/NA for ", pheno,
            ". This script assumes binary phenotypes.")
  }

  # 6.5 Merge with agesex
  pheno_dat <- merge(pheno_dat, agesex_data, by = "IID")
  if (nrow(pheno_dat) == 0) {
    message("  ✗ No overlap between phenotype and agesex data for ", pheno, ". Skipping.")
    next
  }

  # 6.6 Merge all data: PGIs + phenotype + PCs
  merged_data <- merge(
    merge(parental_proband_pgi, pheno_dat, by.x = "true_id", by.y = "IID"),
    pcs,
    by = "true_id"
  )

  if (nrow(merged_data) == 0) {
    message("  ✗ No rows after merging PGIs, phenotype and PCs for ", pheno, ". Skipping.")
    next
  }

  # 6.7 Grouping factor for random intercept
  merged_data$true_id <- factor(merged_data$true_id)

  # 6.8 Rescale continuous covariates & PCs
  cont_vars <- c("age", "age2", "BYEAR", "sexage", "sexage2", paste0("PC", 1:10))
  to_scale  <- intersect(cont_vars, colnames(merged_data))
  if (length(to_scale) > 0) {
    merged_data[to_scale] <- lapply(
      merged_data[to_scale],
      function(x) as.numeric(scale(x))
    )
  }

  # Ensure phenotype is 0/1
  merged_data$phenotype <- as.numeric(merged_data$phenotype)

  # 6.9 Build GLMM formulas (logistic)
  pc_terms  <- paste(paste0("PC", 1:10), collapse = " + ")
  add_terms <- "sex + age + age2 + sexage + sexage2"

  base_rhs <- paste("proband", pc_terms, add_terms, sep = " + ")
  formula_pop_effect    <- as.formula(
    paste("phenotype ~", base_rhs, "+ (1 | true_id)")
  )
  formula_within_family <- as.formula(
    paste("phenotype ~ proband + parental +", pc_terms, "+", add_terms, "+ (1 | true_id)")
  )

  message("  Population model:    ", deparse(formula_pop_effect))
  message("  Within-family model: ", deparse(formula_within_family))

  # 6.10 Fit GLMMs (population + within-family)
  ctrl <- lme4::glmerControl(
    optimizer       = "nloptwrap",
    optCtrl         = list(maxfun = 1e6),
    check.conv.grad = "ignore",
    check.conv.hess = "ignore"
  )

  # Population effect
  message("  → Fitting PopEffect GLMM for ", pheno)
  pop_model <- tryCatch(
    lme4::glmer(
      formula_pop_effect,
      data    = merged_data,
      family  = binomial(),
      control = ctrl
    ),
    error = function(e) {
      message("    ✗ ERROR in PopEffect glmer(): ", e$message)
      NULL
    }
  )

  # Retry with Laplace if needed
  if (!is.null(pop_model)) {
    conv_msg <- pop_model@optinfo$conv$lme4$messages
  } else {
    conv_msg <- NULL
  }

  if (is.null(pop_model) ||
      (!is.null(conv_msg) &&
       grepl("Downdated VtV is not positive definite",
             paste(conv_msg, collapse = " ")))) {
    message("    → Retrying PopEffect with Laplace (nAGQ = 0)")
    pop_model <- tryCatch(
      lme4::glmer(
        formula_pop_effect,
        data    = merged_data,
        family  = binomial(),
        nAGQ    = 0,
        control = ctrl
      ),
      error = function(e) {
        message("    ✗ ERROR in PopEffect glmer() with nAGQ = 0: ", e$message)
        NULL
      }
    )
  }

  if (is.null(pop_model)) {
    message("    ✗ PopEffect GLMM failed for ", pheno, "; skipping phenotype.")
    next
  }

  if (lme4::isSingular(pop_model)) {
    message("    ⚠ PopEffect GLMM singular fit; extracting fixed effects anyway.")
  }

  fe_pop <- summary(pop_model)$coefficients
  if (!"proband" %in% rownames(fe_pop)) {
    message("    ✗ 'proband' coefficient not found in PopEffect model for ", pheno, ". Skipping.")
    next
  }

  coef_pop <- fe_pop["proband", "Estimate"]
  se_pop   <- fe_pop["proband", "Std. Error"]
  or_pop   <- exp(coef_pop)

  message(sprintf("    • PopEffect: log-OR = %.4f ± %.4f; OR = %.3f",
                  coef_pop, se_pop, or_pop))

  # Within-family effect
  message("  → Fitting WithinFamily GLMM for ", pheno)
  fam_model <- tryCatch(
    lme4::glmer(
      formula_within_family,
      data    = merged_data,
      family  = binomial(),
      control = ctrl
    ),
    error = function(e) {
      message("    ✗ ERROR in WithinFamily glmer(): ", e$message)
      NULL
    }
  )

  if (!is.null(fam_model)) {
    conv_msg_fam <- fam_model@optinfo$conv$lme4$messages
  } else {
    conv_msg_fam <- NULL
  }

  if (is.null(fam_model) ||
      (!is.null(conv_msg_fam) &&
       grepl("Downdated VtV is not positive definite",
             paste(conv_msg_fam, collapse = " ")))) {
    message("    → Retrying WithinFamily with Laplace (nAGQ = 0)")
    fam_model <- tryCatch(
      lme4::glmer(
        formula_within_family,
        data    = merged_data,
        family  = binomial(),
        nAGQ    = 0,
        control = ctrl
      ),
      error = function(e) {
        message("    ✗ ERROR in WithinFamily glmer() with nAGQ = 0: ", e$message)
        NULL
      }
    )
  }

  if (is.null(fam_model)) {
    message("    ✗ WithinFamily GLMM failed for ", pheno, "; skipping phenotype.")
    next
  }

  if (lme4::isSingular(fam_model)) {
    message("    ⚠ WithinFamily GLMM singular fit; extracting fixed effects anyway.")
  }

  fe_fam <- summary(fam_model)$coefficients
  if (!all(c("proband", "parental") %in% rownames(fe_fam))) {
    message("    ✗ 'proband' and/or 'parental' coefficients missing in WithinFamily model for ",
            pheno, ". Skipping.")
    next
  }

  coef_fam <- fe_fam["proband",  "Estimate"]
  se_fam   <- fe_fam["proband",  "Std. Error"]
  or_fam   <- exp(coef_fam)

  coef_par <- fe_fam["parental", "Estimate"]
  se_par   <- fe_fam["parental", "Std. Error"]
  or_par   <- exp(coef_par)

  message(sprintf("    • WithinFamily proband: log-OR = %.4f ± %.4f; OR = %.3f",
                  coef_fam, se_fam, or_fam))
  message(sprintf("    • WithinFamily parental: log-OR = %.4f ± %.4f; OR = %.3f",
                  coef_par, se_par, or_par))

  # Sample size used in within-family model
  n_used <- sum(!is.na(lme4::getME(fam_model, "y")))

  # 6.11 Store results for this phenotype
  results_list[[pheno]] <- data.frame(
    Phenotype           = pheno,
    Coef_PopEffect      = coef_pop,
    SE_PopEffect        = se_pop,
    OddsRatio_PopEffect = or_pop,
    Coef_Proband        = coef_fam,
    SE_Proband          = se_fam,
    OddsRatio_Proband   = or_fam,
    Coef_Parental       = coef_par,
    SE_Parental         = se_par,
    OddsRatio_Parental  = or_par,
    N                   = n_used,
    stringsAsFactors    = FALSE
  )

  message("  ✓ Completed GLMM analysis for ", pheno)
}

# ── 7. Write results ─────────────────────────────────────────────────────────

if (length(results_list) == 0) {
  warning("No phenotypes were successfully processed. No output file written.")
} else {
  results_df <- dplyr::bind_rows(results_list)
  readr::write_csv(results_df, results_file)
  message("All results saved to: ", results_file)
}