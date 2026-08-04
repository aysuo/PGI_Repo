#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.6.2_FigS6_IncR2_ratio_byDemography_UKB.R
#
# Purpose:
#   Recreate Figure S6:
#   "Ratio of PGI predictive accuracy (incremental R²) in subgroups of UKB for
#    two traits, comparing models without parental PGI controls (blue circles)
#    to models with parental PGI controls (orange triangles). Panel A shows the
#    female / male incremental R² ratio for diastolic blood pressure; Panel B
#    shows the younger / older incremental R² ratio for body mass index.
#    Error bars are 95% confidence intervals obtained by the percentile method
#    from 1,000 non-parametric bootstrap resamples of families within each
#    subgroup. The dashed horizontal line at 1 indicates equal predictive
#    accuracy across subgroups."
#
# Overview:
#   1) For BPdia and BMI in UKB1, fit two mixed models within demographic
#      subgroups:
#        • "Population Association" model: proband PGI only (no parental PGIs)
#        • "Causal Effect" model: proband + parental PGIs
#
#   2) Compute incremental R² for the proband PGI within each subgroup:
#        BPdia: Male vs Female
#        BMI  : Younger vs Older (median-based age split)
#
#   3) Use 1,000 non-parametric bootstrap resamples of families (FID) to
#      estimate the ratio of R² (group2/group1) and its 95% percentile CI for
#      each model type.
#
#   4) Generate a ratio-only plot (Stratified_R2_Ratio_Only_Plot.png) with:
#        • X-axis: model type (without vs with parental PGIs)
#        • Y-axis: R² ratio, with dashed line at 1
#        • Blue circles: without parental PGIs
#        • Orange triangles: with parental PGIs
#
# Required environment variables (set via paths_PGIrepo_withinFam):
#   PGI_RepoV2            Base project dir  
#   UKB1_PHENO_LIST       Text file with UKB1 phenotypes (one per line); intersected with c("BPdia","BMI")
#   UKB1_COVAR_FILE       Space-delimited covariate file with IID, Sex, PC1–PC20
#   UKB1_PARENTAL_PGS_DIR Dir with PGS_UKB1_parental_<PHENO>-single_SBayesR.pgs.txt
#   UKB1_BMI_PHENO_FILE   BMI noresid phenotype file with IID, BMI (and optionally FID)
#   UKB1_BPdia_PHENO_FILE BPdia noresid phenotype file with IID, BPdia (and optionally FID)
#   UKB1_BPdia_AGE_FILE   File with IID and a single BP-age column (first non-IID column used as age)
#
# Optional:
#   UKB1_R2_RATIO_DIR     Output dir for bootstrap CSVs (default: <PGI_RepoV2>/processed/sumStats/within_family_related)
#   UKB1_FIG_DIR          Output dir for figures       (default: <PGI_RepoV2>/processed/7_Figures/direct_indirect_PGIs)
#
# Output:
#   • Bootstrap_R2_Ratio_Results_BPdia.csv
#   • Bootstrap_R2_Ratio_Results_BMI.csv
#   • Stratified_R2_Ratio_Only_Plot.png
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0) Resolve paths from environment
# ------------------------------------------------------------------------------

PGI_RepoV2 <- Sys.getenv("PGI_RepoV2")
if (PGI_RepoV2 == "") {
  stop("PGI_RepoV2 is not defined. Please set it in your environment (e.g., via paths_PGIrepo_withinFam).")
}

phenolist_file <- Sys.getenv("UKB1_PHENO_LIST")
covar_file     <- Sys.getenv("UKB1_COVAR_FILE")
pgs_dir        <- Sys.getenv("UKB1_PARENTAL_PGS_DIR")
bmi_pheno_file   <- Sys.getenv("UKB1_BMI_PHENO_FILE")
bpdia_pheno_file <- Sys.getenv("UKB1_BPdia_PHENO_FILE")
bpdia_age_file   <- Sys.getenv("UKB1_BPdia_AGE_FILE")

if (phenolist_file == "" || !file.exists(phenolist_file)) {
  stop("UKB1_PHENO_LIST is not set or file does not exist.")
}
if (covar_file == "" || !file.exists(covar_file)) {
  stop("UKB1_COVAR_FILE is not set or file does not exist.")
}
if (pgs_dir == "" || !dir.exists(pgs_dir)) {
  stop("UKB1_PARENTAL_PGS_DIR is not set or directory does not exist.")
}
if (bmi_pheno_file == "" || !file.exists(bmi_pheno_file)) {
  stop("UKB1_BMI_PHENO_FILE is not set or file does not exist.")
}
if (bpdia_pheno_file == "" || !file.exists(bpdia_pheno_file)) {
  stop("UKB1_BPdia_PHENO_FILE is not set or file does not exist.")
}
if (bpdia_age_file == "" || !file.exists(bpdia_age_file)) {
  stop("UKB1_BPdia_AGE_FILE is not set or file does not exist.")
}

ratio_dir_env <- Sys.getenv("UKB1_R2_RATIO_DIR")
if (ratio_dir_env == "") {
  ratio_dir <- file.path(
    PGI_RepoV2,
    "processed", "sumStats", "within_family_related"
  )
} else {
  ratio_dir <- ratio_dir_env
}
if (!dir.exists(ratio_dir)) {
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
}

fig_dir_env <- Sys.getenv("UKB1_FIG_DIR")
if (fig_dir_env == "") {
  fig_dir <- file.path(
    PGI_RepoV2,
    "processed", "7_Figures", "direct_indirect_PGIs"
  )
} else {
  fig_dir <- fig_dir_env
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
}

output_file_ratio_only <- file.path(
  fig_dir,
  "Stratified_R2_Ratio_Only_Plot.png"
)

# ------------------------------------------------------------------------------
# 1) Libraries
# ------------------------------------------------------------------------------

personal_lib_path <- file.path(PGI_RepoV2, "doc", "R", "library")
.libPaths(c(personal_lib_path, .libPaths()))

pkgs <- c("dplyr", "readr", "lme4", "purrr", "tibble", "magrittr", "ggplot2")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}

# ------------------------------------------------------------------------------
# 2) Helper functions
# ------------------------------------------------------------------------------

# Standardization helper
standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Fit LMM/GLMM and compute incremental R² for proband PGI
analyze_effects_stratified <- function(data, formula, family_type = "gaussian") {
  if (nrow(data) < 20) {
    return(NA_real_)
  }

  model <- tryCatch({
    if (family_type == "binomial") {
      glmer(
        formula,
        data    = data,
        family  = binomial(),
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl   = list(maxfun = 2e5))
      )
    } else {
      lmer(
        formula,
        data    = data,
        control = lmerControl(optimizer = "bobyqa",
                              optCtrl   = list(maxfun = 2e5))
      )
    }
  }, error = function(e) NULL)

  if (is.null(model)) return(NA_real_)

  fe <- summary(model)$coefficients
  if (!"proband" %in% rownames(fe)) return(NA_real_)

  beta <- fe["proband", "Estimate"]
  se   <- fe["proband", "Std. Error"]
  if (is.na(beta) || is.na(se) || se == 0) return(NA_real_)

  t_stat <- beta / se
  n_obs  <- nobs(model)
  k_fix  <- length(fixef(model))
  df_res <- n_obs - k_fix - 1
  if (df_res <= 0) return(NA_real_)

  t2  <- t_stat^2
  r2  <- t2 / (t2 + df_res)
  return(r2)
}

# One bootstrap iteration: resample families, compute stratified R² per model
run_bootstrap_iteration_stratified <- function(original_data,
                                               family_ids,
                                               formulas,
                                               strat_var,
                                               group1,
                                               group2) {
  boot_fids <- sample(family_ids, size = length(family_ids), replace = TRUE)
  boot_data <- tibble(FID = boot_fids) %>%
    dplyr::left_join(original_data, by = "FID", relationship = "many-to-many")

  data_g1 <- boot_data %>% dplyr::filter(.data[[strat_var]] == group1)
  data_g2 <- boot_data %>% dplyr::filter(.data[[strat_var]] == group2)

  pop_r2_g1    <- analyze_effects_stratified(data_g1, formulas$pop, formulas$family)
  pop_r2_g2    <- analyze_effects_stratified(data_g2, formulas$pop, formulas$family)
  direct_r2_g1 <- analyze_effects_stratified(data_g1, formulas$direct, formulas$family)
  direct_r2_g2 <- analyze_effects_stratified(data_g2, formulas$direct, formulas$family)

  if (any(is.na(c(pop_r2_g1, pop_r2_g2, direct_r2_g1, direct_r2_g2))) ||
      pop_r2_g1 <= 0 || direct_r2_g1 <= 0) {
    return(tibble(
      pop_r2_g1    = NA_real_,
      pop_r2_g2    = NA_real_,
      direct_r2_g1 = NA_real_,
      direct_r2_g2 = NA_real_
    ))
  }

  tibble(
    pop_r2_g1    = pop_r2_g1,
    pop_r2_g2    = pop_r2_g2,
    direct_r2_g1 = direct_r2_g1,
    direct_r2_g2 = direct_r2_g2
  )
}

# ------------------------------------------------------------------------------
# 3) Bootstrap R² ratios for BPdia and BMI
# ------------------------------------------------------------------------------

# Phenotype list (filter to BPdia and BMI)
phenotypes_all <- read.table(
  phenolist_file,
  header = FALSE,
  stringsAsFactors = FALSE
)$V1
phenotypes <- intersect(phenotypes_all, c("BPdia", "BMI"))

if (length(phenotypes) == 0) {
  stop("Neither 'BPdia' nor 'BMI' found in UKB1_PHENO_LIST. Check your phenotype list.")
}

# Load covariates (IID, Sex, PCs)
covs <- readr::read_delim(
  covar_file,
  delim = " ",
  show_col_types = FALSE
) %>%
  dplyr::select(IID, Sex, dplyr::starts_with("PC"))

# Bootstrap settings
B                <- 1000  # iterations
MIN_BOOT_SUCCESS <- 100   # minimum number of successful iterations

# Container to hold ratio results for plotting
ratio_results_all <- list()

for (pheno in phenotypes) {
  message("\n", strrep("─", 60))
  message("  Stratified incremental R² bootstrap for phenotype: ", pheno)
  message(strrep("─", 60))

  # 3.1 PGI file
  pgi_file <- file.path(
    pgs_dir,
    paste0("PGS_UKB1_parental_", pheno, "-single_SBayesR.pgs.txt")
  )
  if (!file.exists(pgi_file)) {
    message("  ▶ Skipping ", pheno, ": PGI file not found at ", pgi_file)
    next
  }

  # 3.2 Phenotype + age handling
  if (pheno == "BMI") {
    phenotype_file <- bmi_pheno_file
  } else if (pheno == "BPdia") {
    phenotype_file <- bpdia_pheno_file
  } else {
    message("  ▶ Skipping ", pheno, ": phenotype not configured in this script.")
    next
  }

  if (!file.exists(phenotype_file)) {
    message("  ▶ Skipping ", pheno, ": phenotype file not found at ", phenotype_file)
    next
  }

  pgi_data   <- read.table(pgi_file, header = TRUE)
  pheno_data <- read.table(phenotype_file, header = TRUE, fill = TRUE)

  # Drop FID from phenotype file (keep FID from PGI file)
  if ("FID" %in% names(pheno_data)) {
    pheno_data <- pheno_data %>% dplyr::select(-FID)
  }

  if (!pheno %in% names(pheno_data)) {
    stop("Phenotype column ", pheno, " not found in phenotype file: ", phenotype_file)
  }
  pheno_data <- pheno_data %>% dplyr::rename(Phenotype = !!pheno)

  dat <- pgi_data %>%
    dplyr::left_join(pheno_data, by = "IID") %>%
    dplyr::left_join(covs,       by = "IID")

  # Age
  if (pheno == "BMI") {
    if (!"age_at_measurement" %in% names(dat)) {
      stop("Expected column 'age_at_measurement' in BMI data or merged data.")
    }
    dat <- dat %>% dplyr::rename(age = age_at_measurement)
  } else if (pheno == "BPdia") {
    age_data <- read.table(bpdia_age_file, header = TRUE)
    if (!"IID" %in% names(age_data)) {
      stop("BPdia age file must contain 'IID'.")
    }
    age_cols <- setdiff(names(age_data), "IID")
    if (length(age_cols) == 0) {
      stop("BPdia age file must contain at least one age column.")
    }
    age_col <- age_cols[1]
    dat <- dat %>%
      dplyr::left_join(age_data, by = "IID") %>%
      dplyr::rename(age = !!age_col)
  }

  # Clean sample and basic filters
  dat <- dat %>%
    dplyr::filter(
      !is.na(Phenotype),
      !is.na(proband),
      !is.na(Sex),
      !is.na(age),
      !is.na(FID)
    ) %>%
    dplyr::mutate(
      Sex = factor(Sex, levels = 1:2, labels = c("Male", "Female"))
    ) %>%
    droplevels()

  if (nrow(dat) < 100) {
    message("  ▶ Skipping ", pheno, ": too few observations after filtering.")
    next
  }

  # 3.3 Stratification variable and formulas
  pcs <- paste0("PC", 1:20, collapse = " + ")

  if (pheno == "BMI") {
    message("  • Stratifying BMI by AgeGroup (Younger vs Older).")
    med_age <- median(dat$age, na.rm = TRUE)
    dat <- dat %>%
      dplyr::mutate(
        AgeGroup = factor(
          ifelse(age > med_age, "Older", "Younger"),
          levels = c("Younger", "Older")
        )
      )
    strat_var <- "AgeGroup"
    group1    <- "Younger"
    group2    <- "Older"

    form_pop <- as.formula(
      paste0("Phenotype ~ proband + age + Sex + ", pcs, " + (1|FID)")
    )
    form_direct <- as.formula(
      paste0("Phenotype ~ proband + age + Sex + paternal + maternal + ",
             pcs, " + (1|FID)")
    )
  } else {  # BPdia
    message("  • Stratifying BPdia by Sex (Male vs Female).")
    strat_var <- "Sex"
    group1    <- "Male"
    group2    <- "Female"

    form_pop <- as.formula(
      paste0("Phenotype ~ proband + age + ", pcs, " + (1|FID)")
    )
    form_direct <- as.formula(
      paste0("Phenotype ~ proband + age + paternal + maternal + ",
             pcs, " + (1|FID)")
    )
  }

  family_type <- "gaussian"  # both traits are quantitative
  formulas <- list(pop = form_pop, direct = form_direct, family = family_type)

  # 3.4 Standardize variables
  message("  • Standardizing Phenotype within strata and covariates globally.")
  dat <- dat %>%
    dplyr::group_by(.data[[strat_var]]) %>%
    dplyr::mutate(Phenotype = standardize(Phenotype)) %>%
    dplyr::ungroup()

  dat <- dat %>%
    dplyr::mutate(
      proband  = standardize(proband),
      paternal = standardize(paternal),
      maternal = standardize(maternal),
      age      = standardize(age)
    )

  pc_cols <- paste0("PC", 1:20)
  pc_cols <- intersect(pc_cols, names(dat))
  if (length(pc_cols) > 0) {
    dat[pc_cols] <- lapply(dat[pc_cols], standardize)
  }

  # 3.5 Observed R² in each group
  message("  • Computing observed incremental R² in each stratum...")
  observed <- list()
  for (g in c(group1, group2)) {
    gdat <- dat %>% dplyr::filter(.data[[strat_var]] == g)
    observed[[g]] <- list(
      pop    = analyze_effects_stratified(gdat, formulas$pop, formulas$family),
      direct = analyze_effects_stratified(gdat, formulas$direct, formulas$family)
    )
  }

  if (any(is.na(c(
    observed[[group1]]$pop,
    observed[[group2]]$pop,
    observed[[group1]]$direct,
    observed[[group2]]$direct
  )))) {
    message("  ▶ Skipping ", pheno, ": NA observed R² in at least one cell.")
    next
  }

  observed_ratio_pop    <- observed[[group2]]$pop    / observed[[group1]]$pop
  observed_ratio_direct <- observed[[group2]]$direct / observed[[group1]]$direct

  # 3.6 Bootstrap
  message("  • Starting bootstrap with B = ", B, " iterations...")
  unique_fids <- unique(dat$FID)

  boot_results <- purrr::map_dfr(
    seq_len(B),
    function(i) {
      if (i %% 100 == 0) {
        message("    Bootstrap iteration ", i, " / ", B)
      }
      run_bootstrap_iteration_stratified(
        original_data = dat,
        family_ids    = unique_fids,
        formulas      = formulas,
        strat_var     = strat_var,
        group1        = group1,
        group2        = group2
      )
    }
  )

  # Keep successful iterations
  boot_valid <- boot_results %>%
    dplyr::filter(
      !is.na(pop_r2_g1),
      !is.na(pop_r2_g2),
      !is.na(direct_r2_g1),
      !is.na(direct_r2_g2),
      pop_r2_g1 > 0,
      direct_r2_g1 > 0
    )

  n_success <- nrow(boot_valid)
  message("  • Bootstrap complete. ", n_success, "/", B, " iterations succeeded.")

  if (n_success < MIN_BOOT_SUCCESS) {
    message("  ▶ Skipping ", pheno, ": too few successful bootstrap iterations.")
    next
  }

  # Ratios per bootstrap sample
  boot_ratio_pop    <- boot_valid$pop_r2_g2    / boot_valid$pop_r2_g1
  boot_ratio_direct <- boot_valid$direct_r2_g2 / boot_valid$direct_r2_g1

  # Percentile CIs and SEs
  ratio_pop_ci    <- quantile(boot_ratio_pop,    c(0.025, 0.975), na.rm = TRUE)
  ratio_direct_ci <- quantile(boot_ratio_direct, c(0.025, 0.975), na.rm = TRUE)

  se_ratio_pop    <- sd(boot_ratio_pop,    na.rm = TRUE)
  se_ratio_direct <- sd(boot_ratio_direct, na.rm = TRUE)

  ratio_df <- tibble::tibble(
    phenotype        = pheno,
    model            = c("pop", "direct"),
    r2_ratio         = c(observed_ratio_pop, observed_ratio_direct),
    bootstrap_se     = c(se_ratio_pop, se_ratio_direct),
    r2_ratio_ci_low  = c(ratio_pop_ci[1], ratio_direct_ci[1]),
    r2_ratio_ci_high = c(ratio_pop_ci[2], ratio_direct_ci[2])
  )

  # Save per-phenotype CSV
  out_file <- file.path(
    ratio_dir,
    paste0("Bootstrap_R2_Ratio_Results_", pheno, ".csv")
  )
  readr::write_csv(ratio_df, out_file)
  message("  ✔ Saved R² ratio results for ", pheno, " to: ", out_file)

  # Store in memory for plotting
  ratio_results_all[[pheno]] <- ratio_df
}

# ------------------------------------------------------------------------------
# 4) Generate ratio-only plot (Stratified_R2_Ratio_Only_Plot.png)
# ------------------------------------------------------------------------------

if (length(ratio_results_all) == 0) {
  stop("No ratio results available for plotting. Check bootstrap step above.")
}

# Combine results and add panel titles & plot labels
ratio_only_data <- dplyr::bind_rows(ratio_results_all)

# Phenotype-specific metadata
ratio_only_data <- ratio_only_data %>%
  dplyr::mutate(
    panel_title = dplyr::case_when(
      phenotype == "BPdia" ~ "A. Diastolic Blood Pressure",
      phenotype == "BMI"   ~ "B. Body Mass Index",
      TRUE                 ~ phenotype
    ),
    model = dplyr::recode(
      model,
      pop    = "Without Parental PGI Controls",
      direct = "With Parental PGI Controls"
    ),
    model = factor(
      model,
      levels = c("Without Parental PGI Controls",
                 "With Parental PGI Controls")
    )
  )

if (!all(c("BPdia", "BMI") %in% unique(ratio_only_data$phenotype))) {
  warning("Ratio results are not available for both BPdia and BMI. The plot will include only available phenotypes.")
}

# Ratio-only plot
ratio_only_plot <- ggplot2::ggplot(
  ratio_only_data,
  ggplot2::aes(
    x     = model,
    y     = r2_ratio,
    color = model,
    shape = model
  )
) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype   = "dashed",
    color      = "grey40"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = r2_ratio_ci_low, ymax = r2_ratio_ci_high),
    width    = 0.2,
    size     = 0.8,
    position = ggplot2::position_dodge(width = 0.6)
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.6),
    size     = 3
  ) +
  ggplot2::facet_wrap(
    ~ panel_title,
    ncol = 1
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0.8, NA),
    breaks = function(lims) sort(unique(c(1, scales::pretty_breaks()(lims))))
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Without Parental PGI Controls" = "#0072B2",  # blue
      "With Parental PGI Controls"    = "#D55E00"   # orange
    )
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "Without Parental PGI Controls" = 16,  # circle
      "With Parental PGI Controls"    = 17   # triangle
    )
  ) +
  ggplot2::labs(
    x     = NULL,
    y     = "Ratio of Incremental R² (95% Bootstrap CI)",
    color = "",
    shape = ""
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    legend.position   = "bottom",
    strip.background  = ggplot2::element_rect(fill = "grey90", color = "grey50"),
    strip.text        = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x       = ggplot2::element_blank(),
    axis.ticks.x      = ggplot2::element_blank(),
    plot.title        = ggplot2::element_blank(),
    plot.subtitle     = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = output_file_ratio_only,
  plot     = ratio_only_plot,
  width    = 6,
  height   = 5,
  dpi      = 300
)

cat("\n✔ Stratified R² ratio-only plot saved to:\n  ",
    output_file_ratio_only, "\n")