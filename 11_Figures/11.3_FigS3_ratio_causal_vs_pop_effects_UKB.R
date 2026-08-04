#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.5.2_FigS3_ratio_causal_vs_pop_effects_UKB.R
# Purpose: Recreate Figure S3 showing the ratio of causal effects to population
#          associations for PGIs among binary (Panel A) and quantitative
#          (Panel B) phenotypes in UKB.
#
# Description:
#   - Causal effects are estimated in samples of first-degree relatives
#     (family-based mixed models).
#   - Population associations are estimated in samples of unrelated individuals
#     (UKB3 partition).
#   - Ratios are defined as:
#       * For binary traits: ratio of log-odds effects
#           Ratio = log(OR_causal) / log(OR_population)
#       * For quantitative traits: ratio of standardized coefficients
#           Ratio = beta_causal / beta_population
#   - Standard errors of the ratios are computed using the delta method.
#   - Panels:
#       * Panel A: Binary traits (odds ratio scale for the underlying effects;
#                  ratio on log-effect scale).
#       * Panel B: Quantitative traits (standardized linear effects).
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
#   - <PGI_RepoV2>/figures/FigS3A_ratio_binary.png
#   - <PGI_RepoV2>/figures/FigS3B_ratio_quantitative.png
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
figS3A_file <- file.path(fig_output_dir, "FigS3A_ratio_binary.png")
figS3B_file <- file.path(fig_output_dir, "FigS3B_ratio_quantitative.png")

## 3. Load meta-analyzed effect estimates --------------------------------------

# Quantitative (used later for Panel B)
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

# Binary (used for Panel A)
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

# Category labels with line breaks
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

## 6. Panel B: Quantitative traits – ratio of causal to population -------------
# (Effects on standardized linear scale)

# ALZ treated as binary; exclude from quantitative panel
quantitative_data <- merged_quant %>%
  dplyr::filter(
    Phenotype %in% continuous_phenotypes,
    Phenotype != "ALZ"
  )

# Ratio and SE via delta method (direct / population)
quant_ratio <- quantitative_data %>%
  dplyr::mutate(
    delta = Coef_Proband,
    psi   = Coefficient,
    se_delta = SE_Proband,
    se_psi   = StandardError
  ) %>%
  dplyr::mutate(
    Ratio    = delta / psi,
    SE_ratio = sqrt(
      (se_delta / psi)^2 +
        (delta * se_psi / (psi^2))^2
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

p_figS3B <- ggplot2::ggplot(
  quant_ratio,
  ggplot2::aes(x = FullPhenoName, y = Ratio)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Ratio - 1.96 * SE_ratio,
                 ymax = Ratio + 1.96 * SE_ratio),
    width     = 0.3,
    linewidth = 0.6
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype   = "dashed"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Causal Effect : Population Association"
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales = "free_x",
    space  = "free_x"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = "grey80", size = 0.8),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.border     = ggplot2::element_rect(colour = "lightgrey", fill = NA, size = 0.8),
    panel.grid.major = ggplot2::element_line(color = "lightgrey", size = 0.4),
    panel.grid.minor = ggplot2::element_line(color = "lightgrey", size = 0.2),
    axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1),
    strip.placement  = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x     = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = figS3B_file,
  plot     = p_figS3B,
  width    = 14,
  height   = 8
)

cat("Figure S3 Panel B (quantitative traits: ratio) saved to:\n", figS3B_file, "\n")

## 7. Panel A: Binary traits – ratio of log-odds effects -----------------------
# (Effects are ORs; ratio is defined on log-OR scale)

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

# Drop ADHD if desired (too imprecisely estimated)
effects_combined <- effects_combined %>%
  dplyr::filter(Phenotype != "ADHD")

# Compute ratios on log-odds scale
binary_ratio <- effects_combined %>%
  dplyr::mutate(
    # ORs to log-ORs
    logOR_direct = log(Coef_Direct),
    logOR_pop    = log(Coef_Pop),
    se_logOR_direct = SE_Direct,  # SEs are on log-OR scale
    se_logOR_pop    = SE_Pop
  ) %>%
  dplyr::mutate(
    Ratio    = logOR_direct / logOR_pop,
    SE_ratio = abs(Ratio) * sqrt(
      (se_logOR_direct / logOR_direct)^2 +
        (se_logOR_pop    / logOR_pop)^2
    )
  ) %>%
  dplyr::mutate(
    FullPhenoName = trimws(FullPhenoName),
    Phenotype     = factor(Phenotype, levels = ordered_phenotypes)
  ) %>%
  dplyr::arrange(Phenotype) %>%
  dplyr::mutate(
    FullPhenoName = factor(
      FullPhenoName,
      levels = unique(FullPhenoName[order(match(Phenotype, ordered_phenotypes))])
    ),
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

p_figS3A <- ggplot2::ggplot(
  binary_ratio,
  ggplot2::aes(x = FullPhenoName, y = Ratio)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Ratio - 1.96 * SE_ratio,
                 ymax = Ratio + 1.96 * SE_ratio),
    width     = 0.3,
    linewidth = 0.6
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype   = "dashed"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Causal Effect : Population Association (log-OR scale)"
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales = "free_x",
    space  = "free_x",
    labeller = ggplot2::as_labeller(category_labels)
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = "grey80", size = 0.8),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.border     = ggplot2::element_rect(colour = "lightgrey", fill = NA, size = 0.8),
    panel.grid.major = ggplot2::element_line(color = "lightgrey", size = 0.4),
    panel.grid.minor = ggplot2::element_line(color = "lightgrey", size = 0.2),
    axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1),
    strip.placement  = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x     = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = figS3A_file,
  plot     = p_figS3A,
  width    = 14,
  height   = 8
)

cat("Figure S3 Panel A (binary traits: ratio) saved to:\n", figS3A_file, "\n")