#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.5.1_Fig3_causal_vs_population_effects_UKB.R
# Purpose: Recreate Figure 3 Panels A and B showing causal effects and
#          population associations of PGIs for binary and quantitative traits in UKB.
#
# Description:
#   - Panel A (binary traits):
#       Odds ratios from logistic models (population associations) and
#       generalized linear mixed models with a family random intercept (causal).
#   - Panel B (quantitative traits):
#       Standardized coefficients from linear models (population associations)
#       and linear mixed models with a family random intercept (causal).
#
#   This script:
#     (1) Loads meta-analyzed causal and population effect estimates.
#     (2) Classifies phenotypes as binary vs. quantitative 
#     (3) Merges in human-readable phenotype names and phenotype categories.
#     (4) Produces two faceted plots:
#           - Fig3A_proband_vs_population_binary.png
#           - Fig3B_proband_vs_population_quantitative.png
#
# Usage:
#   - Before running, source the path helper:
#         source("paths_PGIrep_withinFam.R")
#
#   - The helper should define at least:
#         PGI_RepoV2        : base directory of the project
#         UKB1_OUTPUT_DIR   : directory with UKB1 causal effect results
#         ukb3_output_dir   : directory with UKB3 population association results
#         ukb1_pheno_dir    : directory with UKB1 phenotype files
#         pheno_mapping_file: file with short ↔ full phenotype names
#
# Output:
#   - <PGI_RepoV2>/figures/Fig3A_proband_vs_population_binary.png
#   - <PGI_RepoV2>/figures/Fig3B_proband_vs_population_quantitative.png
#
# Reference:
#   Alemu et al. (2025).
#   "An Updated Polygenic Index Repository: Expanded Phenotypes,
#    New Cohorts, and Improved Causal Inference."
# ------------------------------------------------------------------------------

## 0. Optional: source path helper ---------------------------------------------
if (file.exists("paths_PGIrep_withinFam.R")) {
  source("paths_PGIrep_withinFam.R")
}

# Ensure PGI_RepoV2 is defined either as R object or environment variable
if (!exists("PGI_RepoV2")) {
  PGI_RepoV2 <- Sys.getenv("PGI_RepoV2")
}
if (is.null(PGI_RepoV2) || PGI_RepoV2 == "") {
  stop("PGI_RepoV2 is not defined. Set it in 'paths_PGIrep_withinFam.R' or as an environment variable.")
}

## 1. Library path and packages -------------------------------------------------

personal_lib_path <- file.path(PGI_RepoV2, "doc/R/library")
.libPaths(c(personal_lib_path, .libPaths()))

required_packages <- c(
  "dplyr",
  "readr",
  "tidyr",
  "ggplot2",
  "forcats",
  "stringr"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  } else {
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}

## 2. Input and output paths ----------------------------------------------------

# Expect these from paths_PGIrep_withinFam.R
if (!exists("UKB1_OUTPUT_DIR") ||
    !exists("ukb3_output_dir") ||
    !exists("ukb1_pheno_dir") ||
    !exists("pheno_mapping_file")) {
  stop("Required path objects (UKB1_OUTPUT_DIR, ukb3_output_dir, ukb1_pheno_dir, pheno_mapping_file) are not defined.\n",
       "Please source 'paths_PGIrep_withinFam.R' before running this script.")
}

# Output directory and files for the figures
fig_output_dir <- file.path(PGI_RepoV2, "figures")
if (!dir.exists(fig_output_dir)) {
  dir.create(fig_output_dir, recursive = TRUE, showWarnings = FALSE)
}
fig3A_file <- file.path(fig_output_dir, "Fig3A_proband_vs_population_binary.png")
fig3B_file <- file.path(fig_output_dir, "Fig3B_proband_vs_population_quantitative.png")

## 3. Load meta-analyzed effect estimates --------------------------------------

# Quantitative (Panel B) – standardized coefficients
popEffect_quant <- readr::read_delim(
  file = file.path(ukb3_output_dir, "popEffect_coefficients.txt"),
  delim = "\t",
  col_types = cols()
)

probandEffect_quant <- readr::read_delim(
  file = file.path(UKB1_OUTPUT_DIR, "compiled_proband_effects.txt"),
  delim = "\t",
  col_types = cols()
)

merged_quant <- dplyr::inner_join(
  popEffect_quant,
  probandEffect_quant,
  by = "Phenotype"
)

# Binary (Panel A) – odds ratios
direct_effects_bin <- readr::read_csv(
  file.path(UKB1_OUTPUT_DIR, "directEffect_UKB1_noresidPheno_coefficients2.csv"),
  show_col_types = FALSE
)

population_effects_bin <- readr::read_delim(
  file = file.path(ukb3_output_dir, "popEffect_noresidPheno_coefficients3.txt"),
  delim = "\t",
  col_types = cols()
)

# Standardize phenotype names in binary population effects (drop numeric suffixes)
population_effects_bin <- population_effects_bin %>%
  dplyr::mutate(
    Phenotype = gsub("[0-9]+$", "", Phenotype)
  )

## 4. Classify phenotypes as binary vs quantitative ----------------------------

phenotypes_all <- unique(c(
  merged_quant$Phenotype,
  direct_effects_bin$Phenotype,
  population_effects_bin$Phenotype
))

binary_phenotypes     <- character(0)
continuous_phenotypes <- character(0)

for (pheno in phenotypes_all) {
  pheno_file <- file.path(ukb1_pheno_dir, paste0(pheno, "_noresid.pheno"))

  if (!file.exists(pheno_file)) {
    message("Phenotype file not found (skipping classification): ", pheno_file)
    next
  }

  lines <- readLines(pheno_file)
  lines <- gsub("[ \t]+", "\t", lines)

  pheno_data <- utils::read.table(
    text = lines,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    na.strings = "NA",
    quote = "",
    dec = ".",
    stringsAsFactors = FALSE
  )

  if (ncol(pheno_data) < 3) {
    message("Unexpected number of columns in phenotype file (skipping): ", pheno_file)
    next
  }

  pheno_values <- stats::na.omit(pheno_data[[3]])
  if (length(pheno_values) == 0) {
    next
  }

  if (all(pheno_values %in% c(0, 1))) {
    binary_phenotypes <- c(binary_phenotypes, pheno)
  } else {
    continuous_phenotypes <- c(continuous_phenotypes, pheno)
  }
}

binary_phenotypes     <- unique(binary_phenotypes)
continuous_phenotypes <- unique(continuous_phenotypes)

## 5. Phenotype name mapping and categories ------------------------------------

# Short ↔ full phenotype names
pheno_mapping <- readr::read_delim(
  file = pheno_mapping_file,
  delim = "\t",
  col_names = c("Phenotype", "FullPhenoName"),
  col_types = cols()
)

# Phenotype categories (short names)
phenotype_categories <- list(
  "Anthropometric"                 = c("BMI", "HEIGHT"),
  "Blood Biomarkers"               = c("BL_HDL", "BL_LDL", "BL_nonHDL", "BL_CHOL", "BL_TRYG",
                                       "BPdia", "BPpulse", "BPsys"),
  "Cognition & Education"          = c("CP", "EA", "ALZnoproxy"),
  "Personality & Wellbeing"        = c("ADVENTURE", "EXTRA", "FAMSAT", "FRIENDSAT", "MORNING",
                                       "NEURO", "OPEN", "ACTIVITY", "RELIGATT", "RISK", "SWB"),
  "Health"                         = c("ALLERGYPOLLEN", "ASTHMA", "ASTECZRHI", "BRCA", "COPD",
                                       "CAD", "HAYFEVER", "IBD", "MIGRAINE", "NEARSIGHTED",
                                       "PRCA", "SELFHEALTH", "T2D"),
  "Fertility & Sexual Development" = c("AFB", "MENARCHE", "AFS", "NEBwomen"),
  "Psychiatric Conditions"         = c("ADHD", "ANOREX", "ASD", "BIPOLAR", "DEP", "INSOMNIA", "SCZ"),
  "Substance Use"                  = c("AUDIT", "ASI", "CANNABIS", "CPD", "DPW", "EVERSMOKE", "SMCESS")
)

# Category mapping
category_map <- utils::stack(phenotype_categories) %>%
  dplyr::rename(Category = ind, Phenotype = values) %>%
  dplyr::distinct()

pheno_mapping <- pheno_mapping %>%
  dplyr::left_join(category_map, by = "Phenotype")

# Category labels with line breaks for facet strips
category_labels <- c(
  "Anthropometric"                 = "Anthropometric",
  "Blood Biomarkers"               = "Blood\nBiomarkers",
  "Cognition & Education"          = "Cognition &\nEducation",
  "Personality & Wellbeing"        = "Personality &\nWellbeing",
  "Health"                         = "Health",
  "Fertility & Sexual Development" = "Fertility &\nSexual Development",
  "Psychiatric Conditions"         = "Psychiatric\nConditions",
  "Substance Use"                  = "Substance Use"
)

# Ordered full names following phenotype_categories ordering
ordered_full_names <- unlist(lapply(phenotype_categories, function(cat_phenos) {
  pheno_mapping %>%
    dplyr::filter(Phenotype %in% cat_phenos) %>%
    dplyr::arrange(match(Phenotype, cat_phenos)) %>%
    dplyr::pull(FullPhenoName)
}))
ordered_full_names <- unique(ordered_full_names)

# Ordered short names
ordered_phenotypes <- unique(unlist(phenotype_categories))

## 6. Panel B: Quantitative traits (standardized effects) ----------------------

# ALZ treated as binary in the manuscript; exclude from quantitative panel
quantitative_data <- merged_quant %>%
  dplyr::filter(
    Phenotype %in% continuous_phenotypes,
    Phenotype != "ALZ"
  )

# Long format for plotting quantitative effects
long_quant <- quantitative_data %>%
  tidyr::pivot_longer(
    cols      = c(Coef_Proband, Coefficient),
    names_to  = "CoefficientType",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    SE = dplyr::case_when(
      CoefficientType == "Coefficient"   ~ StandardError,
      CoefficientType == "Coef_Proband"  ~ SE_Proband
    ),
    EffectType = dplyr::if_else(
      CoefficientType == "Coef_Proband",
      "Causal Effect",
      "Population Association"
    )
  ) %>%
  dplyr::left_join(pheno_mapping, by = "Phenotype") %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      Phenotype %in% phenotype_categories$Anthropometric                 ~ "Anthropometric",
      Phenotype %in% phenotype_categories$`Blood Biomarkers`             ~ "Blood Biomarkers",
      Phenotype %in% phenotype_categories$`Cognition & Education`        ~ "Cognition & Education",
      Phenotype %in% phenotype_categories$`Personality & Wellbeing`      ~ "Personality & Wellbeing",
      Phenotype %in% phenotype_categories$Health                         ~ "Health",
      Phenotype %in% phenotype_categories$`Fertility & Sexual Development` ~ "Fertility & Sexual Development",
      Phenotype %in% phenotype_categories$`Psychiatric Conditions`       ~ "Psychiatric Conditions",
      Phenotype %in% phenotype_categories$`Substance Use`                ~ "Substance Use",
      TRUE                                                                ~ NA_character_
    ),
    FullPhenoName = factor(FullPhenoName, levels = ordered_full_names),
    Category      = factor(
      Category,
      levels = names(category_labels),
      labels = category_labels
    )
  ) %>%
  dplyr::filter(!is.na(Category), !is.na(FullPhenoName))

dodge_width <- ggplot2::position_dodge(width = 0.5)

p_fig3B <- ggplot2::ggplot(
  long_quant,
  ggplot2::aes(x = FullPhenoName, y = Value, color = EffectType, group = EffectType)
) +
  ggplot2::geom_point(position = dodge_width, size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Value - SE, ymax = Value + SE),
    position = dodge_width,
    width    = 0.4,
    linewidth = 0.6
  ) +
  ggplot2::scale_color_manual(
    values = c("Causal Effect" = "darkblue", "Population Association" = "orange")
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Standardized Effect Size"
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales  = "free_x",
    space   = "free_x"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1),
    axis.title.y     = ggplot2::element_text(face = "plain"),
    legend.position  = "none",
    strip.placement  = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x     = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = fig3B_file,
  plot     = p_fig3B,
  width    = 14,
  height   = 8
)

cat("Figure 3 Panel B (quantitative traits) saved to:\n", fig3B_file, "\n")

## 7. Panel A: Binary traits (odds ratios) -------------------------------------

# Filter and merge datasets for binary traits
binary_direct_effects <- direct_effects_bin %>%
  dplyr::filter(Phenotype %in% binary_phenotypes) %>%
  dplyr::select(
    Phenotype,
    Coef_Direct = Coef_Proband,
    SE_Direct   = Boot_SE_ProbandEffect
  )

binary_population_effects <- population_effects_bin %>%
  dplyr::filter(Phenotype %in% binary_phenotypes) %>%
  dplyr::select(
    Phenotype,
    Coef_Pop = Coef_Proband,
    SE_Pop   = SE_Proband
  )

effects_combined <- binary_direct_effects %>%
  dplyr::inner_join(binary_population_effects, by = "Phenotype") %>%
  dplyr::left_join(pheno_mapping, by = "Phenotype")

# Long format for binary traits, matching coefficients and SEs by effect type
long_bin <- effects_combined %>%
  tidyr::pivot_longer(
    cols      = c(Coef_Direct, Coef_Pop),
    names_to  = "EffectType",
    values_to = "Coefficient"
  ) %>%
  tidyr::pivot_longer(
    cols      = c(SE_Direct, SE_Pop),
    names_to  = "SEType",
    values_to = "SE"
  ) %>%
  dplyr::mutate(
    EffectType = dplyr::recode(
      EffectType,
      "Coef_Direct" = "Causal Effect",
      "Coef_Pop"    = "Population Association"
    ),
    SEType = dplyr::recode(
      SEType,
      "SE_Direct" = "Causal Effect",
      "SE_Pop"    = "Population Association"
    )
  ) %>%
  dplyr::filter(EffectType == SEType) %>%
  dplyr::select(-SEType)

# Optional: drop ADHD (too imprecisely estimated)
long_bin <- long_bin %>%
  dplyr::filter(Phenotype != "ADHD")

# Order phenotypes and full names
long_bin <- long_bin %>%
  dplyr::mutate(
    Phenotype = factor(Phenotype, levels = ordered_phenotypes)
  ) %>%
  dplyr::arrange(Phenotype) %>%
  dplyr::mutate(
    FullPhenoName = factor(
      FullPhenoName,
      levels = unique(FullPhenoName[order(match(Phenotype, ordered_phenotypes))])
    )
  )

# Category assignment
long_bin <- long_bin %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      Phenotype %in% phenotype_categories$Anthropometric                 ~ "Anthropometric",
      Phenotype %in% phenotype_categories$`Blood Biomarkers`             ~ "Blood Biomarkers",
      Phenotype %in% phenotype_categories$`Cognition & Education`        ~ "Cognition & Education",
      Phenotype %in% phenotype_categories$`Personality & Wellbeing`      ~ "Personality & Wellbeing",
      Phenotype %in% phenotype_categories$Health                         ~ "Health",
      Phenotype %in% phenotype_categories$`Fertility & Sexual Development` ~ "Fertility & Sexual Development",
      Phenotype %in% phenotype_categories$`Psychiatric Conditions`       ~ "Psychiatric Conditions",
      Phenotype %in% phenotype_categories$`Substance Use`                ~ "Substance Use",
      TRUE                                                                ~ NA_character_
    ),
    Category = factor(
      Category,
      levels = names(category_labels),
      labels = category_labels
    )
  ) %>%
  dplyr::filter(!is.na(Category), !is.na(FullPhenoName))

# 95% CIs on odds ratio scale (assuming Coefficient is odds ratio and SE on log scale)
long_bin <- long_bin %>%
  dplyr::mutate(
    Lower_CI_log = log(Coefficient) - 1.96 * SE,
    Upper_CI_log = log(Coefficient) + 1.96 * SE,
    Lower_CI     = exp(Lower_CI_log),
    Upper_CI     = exp(Upper_CI_log)
  )

p_fig3A <- ggplot2::ggplot(
  long_bin,
  ggplot2::aes(x = FullPhenoName, y = Coefficient, color = EffectType, group = EffectType)
) +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5), size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Lower_CI, ymax = Upper_CI),
    position = ggplot2::position_dodge(0.5),
    width    = 0.4,
    linewidth = 0.6
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Causal Effect"          = "darkblue",
      "Population Association" = "orange"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Odds Ratio"
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales   = "free_x",
    space    = "free_x",
    labeller = ggplot2::as_labeller(category_labels)
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.background      = ggplot2::element_rect(fill = "white", color = "white"),
    panel.background     = ggplot2::element_rect(fill = "white", color = "white"),
    panel.border         = ggplot2::element_rect(colour = "#848383", fill = NA, size = 0.8),
    panel.grid.major.x   = ggplot2::element_line(color = "lightgrey", size = 0.3),
    panel.grid.major.y   = ggplot2::element_line(color = "lightgrey", size = 0.3),
    panel.grid.minor     = ggplot2::element_blank(),
    strip.placement      = "outside",
    strip.background     = ggplot2::element_blank(),
    strip.text.x         = ggplot2::element_text(size = 11, face = "bold",
                                                 margin = ggplot2::margin(t = 5, b = 5)),
    axis.text.x          = ggplot2::element_text(angle = 60, hjust = 1),
    axis.title.y         = ggplot2::element_text(face = "plain"),
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction     = "vertical",
    legend.title         = ggplot2::element_blank(),
    axis.line            = ggplot2::element_line(color = "black")
  )

ggplot2::ggsave(
  filename = fig3A_file,
  plot     = p_fig3A,
  width    = 14,
  height   = 8,
  limitsize = FALSE
)

cat("Figure 3 Panel A (binary traits) saved to:\n", fig3A_file, "\n")
