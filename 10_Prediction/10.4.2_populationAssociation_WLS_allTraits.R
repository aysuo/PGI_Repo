#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.4.2_populationAssociation_WLS_allTraits.R
# Purpose: Estimate the population association of polygenic indices (PGIs) for
#          all traits (binary and quantitative) in the WLS unrelated sample.
#
# Description:
#   - Constructs an *unrelated* WLS sample by excluding individuals with
#     parental PGIs (identified via HEIGHT PGIs).
#   - For each phenotype in the WLS list:
#       * Reads the phenotype file (binary: noresid; quantitative: resid).
#       * Merges with standard PGIs (restricted-SNPs), agesex covariates, and PCs.
#       * Standardizes PGIs and quantitative outcomes.
#       * Fits:
#           - Binary traits: logistic regression
#               phenotype ~ proband + PCs + sex + age + age2 + sexage + sexage2
#           - Quantitative traits: linear regression (standardized outcome)
#               phenotype ~ proband + PCs + sex + age + age2 + sexage + sexage2
#   - Outputs coefficient, SE, odds ratio (for binary traits), and N.
#
# Usage:
#   1) Ensure the following environment variables are set, e.g. via
#      `paths_PGIrepo_withinFam` (bash):
#
#         PGI_RepoV2            : base project directory
#         WLS_PHENO_LIST        : phenotype list file (one phenotype per line)
#         WLS_AGESEX_FILE       : agesex covariate file
#         WLS_PC_FILE           : principal components file (WLS_PCs.eigenvec)
#         WLS_PGI_STANDARD_DIR  : directory with standard PGIs (restricted SNPs)
#         WLS_PGI_PARENTAL_DIR  : directory with parental PGIs
#         WLS_PHENO_NORESD_DIR  : directory with binary (noresid) phenotypes
#         WLS_PHENO_RESD_DIR    : directory with quantitative (resid) phenotypes
#         WLS_MISC_DIR          : directory for misc outputs (unrelated sample)
#         WLS_OUTPUT_DIR        : directory for population association results
#
#   2) Run:
#        Rscript 10.4.2_populationAssociation_WLS_allTraits.R
#
# Output:
#   - <WLS_MISC_DIR>/unrelated_sample_WLS.txt
#       (one true_id per line; unrelated sample)
#   - <WLS_OUTPUT_DIR>/popEffect_WLS_allTraits_coefficients.csv
#       Columns:
#         Phenotype, Coef_PopEffect, SE_PopEffect, OddsRatio_PopEffect, N
#
# Reference:
#   Alemu et al. (2025). "An Updated Polygenic Index Repository:
#   Expanded Phenotypes, New Cohorts, and Improved Causal Inference."
# ------------------------------------------------------------------------------

# ── 1. Setup ──────────────────────────────────────────────────────────────────

PGI_RepoV2 <- Sys.getenv("PGI_RepoV2")
if (PGI_RepoV2 == "") {
  stop("PGI_RepoV2 is not defined. Please set it in your environment before running.")
}

# WLS-specific paths (from bash env)
phenotype_list_file   <- Sys.getenv("WLS_PHENO_LIST")
agesex_file           <- Sys.getenv("WLS_AGESEX_FILE")
pc_file               <- Sys.getenv("WLS_PC_FILE")
pgi_standard_dir      <- Sys.getenv("WLS_PGI_STANDARD_DIR")
pgi_parental_dir      <- Sys.getenv("WLS_PGI_PARENTAL_DIR")
pheno_n_resid_dir     <- Sys.getenv("WLS_PHENO_NORESD_DIR")
pheno_resid_dir       <- Sys.getenv("WLS_PHENO_RESD_DIR")
wls_misc_dir          <- Sys.getenv("WLS_MISC_DIR")
wls_output_dir        <- Sys.getenv("WLS_OUTPUT_DIR")

required_vars <- c(
  phenotype_list_file,
  agesex_file,
  pc_file,
  pgi_standard_dir,
  pgi_parental_dir,
  pheno_n_resid_dir,
  pheno_resid_dir,
  wls_misc_dir,
  wls_output_dir
)

if (any(required_vars == "")) {
  stop(
    "One or more required WLS path environment variables are not set.\n",
    "Required: WLS_PHENO_LIST, WLS_AGESEX_FILE, WLS_PC_FILE,\n",
    "          WLS_PGI_STANDARD_DIR, WLS_PGI_PARENTAL_DIR,\n",
    "          WLS_PHENO_NORESD_DIR, WLS_PHENO_RESD_DIR,\n",
    "          WLS_MISC_DIR, WLS_OUTPUT_DIR.\n"
  )
}

if (!dir.exists(wls_misc_dir)) {
  dir.create(wls_misc_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(wls_output_dir)) {
  dir.create(wls_output_dir, recursive = TRUE, showWarnings = FALSE)
}

unrelated_file <- file.path(wls_misc_dir, "unrelated_sample_WLS.txt")
results_file   <- file.path(wls_output_dir, "popEffect_WLS_allTraits_coefficients.csv")

# Set library path for reproducible installs
personal_lib_path <- file.path(PGI_RepoV2, "doc/R/library")
.libPaths(c(personal_lib_path, .libPaths()))

# ── 2. Load packages ─────────────────────────────────────────────────────────

packages <- c("ggplot2", "dplyr", "lmtest", "boot", "readr", "tidyverse", "broom")

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}
for (pkg in packages) install_and_load(pkg)

# ── 3. Identify unrelated sample (based on HEIGHT PGIs) ──────────────────────

message("Identifying unrelated sample using HEIGHT PGIs...")

height_pgi_all_file <- file.path(
  pgi_standard_dir,
  "PGS_WLS_HEIGHT-single_SBayesR.txt"
)
height_pgi_parental_file <- file.path(
  pgi_parental_dir,
  "PGS_WLS_parental_HEIGHT-single_SBayesR.pgs.txt"
)

if (!file.exists(height_pgi_all_file)) {
  stop("HEIGHT standard PGI file not found: ", height_pgi_all_file)
}
if (!file.exists(height_pgi_parental_file)) {
  stop("HEIGHT parental PGI file not found: ", height_pgi_parental_file)
}

all_iids <- read.table(height_pgi_all_file, header = TRUE, sep = "\t")
related_iids <- read.table(height_pgi_parental_file, header = TRUE, sep = " ")

all_iids <- all_iids %>%
  dplyr::mutate(true_id = sub("_.*", "", IID))

related_iids <- related_iids %>%
  dplyr::mutate(true_id = sub("_.*", "", IID))

unrelated_iids <- all_iids %>%
  dplyr::filter(!true_id %in% related_iids$true_id) %>%
  dplyr::select(true_id)

message("Number of unrelated IIDs: ", nrow(unrelated_iids))
utils::head(unrelated_iids)

readr::write_delim(
  unrelated_iids,
  unrelated_file,
  col_names = FALSE
)
message("Unrelated sample saved to: ", unrelated_file)

# ── 4. Read phenotypes list, covariates, PCs ─────────────────────────────────

phenotypes <- readr::read_lines(phenotype_list_file)
phenotypes <- phenotypes[phenotypes != ""]
message("Number of phenotypes in WLS list: ", length(phenotypes))

# Agesex covariates
agesex_data <- read.table(
  agesex_file,
  header    = TRUE,
  fill      = TRUE,
  na.strings = ""
)

agesex_data$IID <- sapply(strsplit(as.character(agesex_data$IID), "_"), `[`, 1)

agesex_data <- agesex_data %>%
  dplyr::group_by(IID) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

agesex_data$sex <- as.numeric(agesex_data$sex == "F")
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

# Unrelated IIDs (as data frame with col "true_id")
unrelated_iids <- readr::read_delim(
  unrelated_file,
  col_names = FALSE,
  show_col_types = FALSE
)
colnames(unrelated_iids) <- c("true_id")

# ── 5. Main loop over phenotypes ─────────────────────────────────────────────

results_list <- list()

for (pheno in phenotypes) {
  message("────────────────────────────────────────────")
  message("Processing phenotype: ", pheno)

  # 5.1 Phenotype file (binary: noresid; quantitative: resid)
  pheno_file_nresid <- file.path(pheno_n_resid_dir, paste0(pheno, ".pheno"))
  pheno_file_resid  <- file.path(pheno_resid_dir, paste0(pheno, ".pheno"))

  if (file.exists(pheno_file_nresid)) {
    phenotype_file <- pheno_file_nresid
  } else if (file.exists(pheno_file_resid)) {
    phenotype_file <- pheno_file_resid
  } else {
    message("  ✗ No phenotype file found for ", pheno, ". Skipping.")
    next
  }

  is_binary <- grepl("noresid", phenotype_file)

  # 5.2 Standard PGI file (restricted SNPs)
  if (pheno == "HARDCAD") {
    # No HARDCAD PGI: use CAD
    standard_pgi_file <- file.path(
      pgi_standard_dir,
      "PGS_WLS_CAD-single_SBayesR.txt"
    )
  } else {
    standard_pgi_file <- file.path(
      pgi_standard_dir,
      paste0("PGS_WLS_", pheno, "-single_SBayesR.txt")
    )
  }

  if (!file.exists(standard_pgi_file)) {
    message("  ✗ Standard PGI file not found for ", pheno, ": ", standard_pgi_file,
            ". Skipping.")
    next
  }

  standard_pgi <- read.table(standard_pgi_file, header = TRUE, sep = "\t")

  # Last column is the PGI; rename to "proband"
  colnames(standard_pgi)[ncol(standard_pgi)] <- "proband"
  standard_pgi$true_id <- sub("_.*", "", as.character(standard_pgi$IID))

  # Restrict to unrelated sample
  standard_pgi <- dplyr::inner_join(
    standard_pgi,
    unrelated_iids,
    by = "true_id"
  )

  if (nrow(standard_pgi) == 0) {
    message("  ✗ No overlap between PGI and unrelated sample for ", pheno, ". Skipping.")
    next
  }

  # 5.3 Phenotype data
  pheno_dat <- read.table(
    phenotype_file,
    header = TRUE,
    sep    = ","
  )
  pheno_dat$true_id <- sub("_.*", "", as.character(pheno_dat$id))

  # Binary phenotypes sometimes use "pheno"
  if (is_binary && "pheno" %in% names(pheno_dat)) {
    names(pheno_dat)[names(pheno_dat) == "pheno"] <- "phenotype"
  }

  if (!"phenotype" %in% names(pheno_dat)) {
    message("  ✗ 'phenotype' column not found for ", pheno, ". Skipping.")
    next
  }

  # 5.4 Merge phenotype with agesex
  pheno_dat <- dplyr::inner_join(
    pheno_dat,
    agesex_data,
    by = c("true_id" = "IID")
  )

  if (nrow(pheno_dat) == 0) {
    message("  ✗ No overlap between phenotype and agesex data for ", pheno, ". Skipping.")
    next
  }

  # 5.5 Merge phenotype + covariates + PGI
  merged_dat <- dplyr::inner_join(
    pheno_dat,
    standard_pgi,
    by = "true_id"
  )

  if (nrow(merged_dat) == 0) {
    message("  ✗ No overlap between phenotype+covariates and PGI for ", pheno,
            ". Skipping.")
    next
  }

  # 5.6 Merge with PCs
  merged_dat <- dplyr::inner_join(
    merged_dat,
    pcs,
    by = "true_id"
  )

  if (nrow(merged_dat) == 0) {
    message("  ✗ No overlap between merged data and PCs for ", pheno, ". Skipping.")
    next
  }

  # 5.7 Standardization and checks
  # Standardize PGI
  merged_dat$proband <- as.numeric(scale(merged_dat$proband))

  if (!is_binary) {
    merged_dat$phenotype <- as.numeric(scale(merged_dat$phenotype))
  } else {
    # Ensure binary-coded
    merged_dat$phenotype <- as.numeric(merged_dat$phenotype)
    if (!all(merged_dat$phenotype %in% c(0, 1, NA))) {
      message("  ⚠ Phenotype not strictly 0/1 for ", pheno,
              " (still treated as binary logit).")
    }
  }

  # 5.8 Build regression formula
  pc_terms  <- intersect(paste0("PC", 1:10), colnames(merged_dat))
  cov_terms <- intersect(
    c("sex", "age", "age2", "sexage", "sexage2"),
    colnames(merged_dat)
  )

  rhs_terms <- c("proband", pc_terms, cov_terms)
  if (length(rhs_terms) == 0) {
    message("  ✗ No RHS terms available for regression for ", pheno, ". Skipping.")
    next
  }

  formula_str <- paste("phenotype ~", paste(rhs_terms, collapse = " + "))
  model_formula <- as.formula(formula_str)

  message("  Model formula: ", formula_str)

  # 5.9 Fit model and extract estimates
  if (is_binary) {
    fit <- glm(
      model_formula,
      data   = merged_dat,
      family = binomial(link = "logit")
    )

    if (!"proband" %in% names(coef(fit))) {
      message("  ✗ 'proband' coefficient missing in glm for ", pheno, ". Skipping.")
      next
    }

    coef_est    <- coef(fit)["proband"]
    vcov_mat    <- stats::vcov(fit)
    se_log_odds <- sqrt(diag(vcov_mat)["proband"])
    odds_ratio  <- exp(coef_est)
    n_used      <- stats::nobs(fit)

    se_pop <- se_log_odds

  } else {
    fit <- lm(
      model_formula,
      data = merged_dat
    )

    if (!"proband" %in% names(coef(fit))) {
      message("  ✗ 'proband' coefficient missing in lm for ", pheno, ". Skipping.")
      next
    }

    coef_est <- coef(fit)["proband"]
    se_coef  <- summary(fit)$coefficients["proband", "Std. Error"]
    odds_ratio <- NA_real_
    n_used <- stats::nobs(fit)

    se_pop <- se_coef
  }

  result_row <- data.frame(
    Phenotype           = pheno,
    Coef_PopEffect      = coef_est,
    SE_PopEffect        = se_pop,
    OddsRatio_PopEffect = odds_ratio,
    N                   = n_used,
    stringsAsFactors    = FALSE
  )

  results_list[[pheno]] <- result_row

  message("  ✓ Completed population association for ", pheno)
}

# ── 6. Save results ──────────────────────────────────────────────────────────

if (length(results_list) == 0) {
  warning("No phenotypes were successfully processed. No output file written.")
} else {
  results_df <- dplyr::bind_rows(results_list)
  readr::write_csv(results_df, results_file)
  message("All population association results saved to: ", results_file)
}