#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.5.3_FigS4_causal_vs_population_effects_WLS.R
# Purpose: Recreate Figure S4 showing causal effects of PGIs versus their
#          population associations in the WLS cohort.
#
# Description:
#   This script visualizes:
#     • Panel A (binary phenotypes):
#         - Causal effects: odds ratios from generalized linear mixed models
#           with a logistic link and a random intercept for family ID, estimated
#           in first-degree relatives (sample with sibling genotype data).
#         - Population associations: odds ratios from standard logistic
#           regression, estimated in a sample of unrelated individuals.
#
#     • Panel B (quantitative phenotypes):
#         - Causal effects: standardized coefficients from linear mixed models
#           (random intercept for family ID) in the relative sample.
#         - Population associations: standardized coefficients from linear
#           regression in the unrelated sample.
#
#   This script:
#     (1) Reads pre-computed causal and population association estimates
#         produced by:
#           - 10.4.1_causalEffect_WLS_BinaryTraits_lme4.R
#           - 10.4.2_populationAssociation_WLS_allTraits.R
#     (2) Merges in phenotype names and categories.
#     (3) Generates faceted plots for binary and quantitative traits:
#           - FigS4A_WLS_binary_causal_vs_population.png
#           - FigS4B_WLS_quantitative_causal_vs_population.png
#
# Usage:
#   - Ensure the following environment variables are set (e.g. via
#     `paths_PGIrepo_withinFam` bash file):
#
#       PGI_RepoV2          : base project directory
#       WLS_OUTPUT_DIR      : directory with WLS effect summary files
#       WLS_FIG_DIR         : directory for WLS figures (optional; defaults to
#                             <PGI_RepoV2>/figures if unset)
#       PHENO_MAPPING_FILE  : file with short ↔ full phenotype names
#
#   - Then run:
#       Rscript 10.5.3_FigS4_causal_vs_population_effects_WLS.R
#
# Input:
#   - <WLS_OUTPUT_DIR>/binaryNonBinary_directPGIeffect_coefficients.txt
#       (causal/direct effects from mixed models in relatives)
#   - <WLS_OUTPUT_DIR>/binaryNonBinary_popPGIeffect_coefficients.txt
#       (population associations in unrelated WLS sample)
#   - <PHENO_MAPPING_FILE> (tab-delimited: Phenotype, FullPhenoName)
#
# Output:
#   - <WLS_FIG_DIR>/FigS4A_WLS_binary_causal_vs_population.png
#   - <WLS_FIG_DIR>/FigS4B_WLS_quantitative_causal_vs_population.png
#
# Reference:
#   Alemu et al. (2025). "An Updated Polygenic Index Repository:
#   Expanded Phenotypes, New Cohorts, and Improved Causal Inference."
# ------------------------------------------------------------------------------

## 1. Setup --------------------------------------------------------------------

PGI_RepoV2 <- Sys.getenv("PGI_RepoV2")
if (PGI_RepoV2 == "") {
  stop("PGI_RepoV2 is not defined. Please set it in your environment before running.")
}

wls_output_dir     <- Sys.getenv("WLS_OUTPUT_DIR")
wls_fig_dir_env    <- Sys.getenv("WLS_FIG_DIR")
pheno_mapping_file <- Sys.getenv("PHENO_MAPPING_FILE")

if (wls_output_dir == "") {
  stop("WLS_OUTPUT_DIR is not defined. Set it in your environment (paths_PGIrepo_withinFam).")
}
if (pheno_mapping_file == "") {
  stop("PHENO_MAPPING_FILE is not defined. Set it in your environment (paths_PGIrepo_withinFam).")
}

# Default figure directory if not explicitly set
if (wls_fig_dir_env == "") {
  wls_fig_dir <- file.path(PGI_RepoV2, "figures")
} else {
  wls_fig_dir <- wls_fig_dir_env
}
if (!dir.exists(wls_fig_dir)) {
  dir.create(wls_fig_dir, recursive = TRUE, showWarnings = FALSE)
}

figS4A_file <- file.path(wls_fig_dir, "FigS4A_WLS_binary_causal_vs_population.png")
figS4B_file <- file.path(wls_fig_dir, "FigS4B_WLS_quantitative_causal_vs_population.png")

# Effect summary files (created by the WLS causal & population scripts)
direct_file <- file.path(
  wls_output_dir,
  "binaryNonBinary_directPGIeffect_coefficients.txt"
)
pop_file <- file.path(
  wls_output_dir,
  "binaryNonBinary_popPGIeffect_coefficients.txt"
)

## 2. Libraries ----------------------------------------------------------------

personal_lib_path <- file.path(PGI_RepoV2, "doc/R/library")
.libPaths(c(personal_lib_path, .libPaths()))

packages <- c("dplyr", "readr", "tidyr", "ggplot2", "forcats", "stringr")

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}
for (pkg in packages) install_and_load(pkg)

## 3. Phenotype mapping and categories ----------------------------------------

# Short ↔ full phenotype names (tab-delimited)
pheno_mapping <- readr::read_delim(
  file      = pheno_mapping_file,
  delim     = "\t",
  col_names = c("Phenotype", "FullPhenoName"),
  col_types = cols()
)

# Phenotype categories (short names)
phenotype_categories <- list(
  "Anthropometric"                 = c("BMI", "HEIGHT"),
  "Blood Biomarkers"               = c("BL_HDL", "BL_LDL", "BL_nonHDL", "BL_CHOL", "BL_TRYG",
                                       "BPdia", "BPpulse", "BPsys"),
  "Cognition & Education"          = c("CP", "EA", "ALZnoproxy", "ALZ"),
  "Personality & Wellbeing"        = c("ADVENTURE", "EXTRA", "FAMSAT", "FRIENDSAT", "MORNING",
                                       "NEURO", "OPEN", "ACTIVITY", "RELIGATT", "RISK", "SWB"),
  "Health"                         = c("ALLERGYPOLLEN", "ASTHMA", "ASTECZRHI", "BRCA", "COPD",
                                       "HARDCAD", "HAYFEVER", "IBD", "MIGRAINE", "NEARSIGHTED",
                                       "PRCA", "SELFHEALTH", "T2D"),
  "Fertility & Sexual Development" = c("AFB", "MENARCHE", "AFS", "NEBwomen"),
  "Psychiatric Conditions"         = c("ADHD", "ANOREX", "ASD", "BIPOLAR", "DEP", "INSOMNIA", "SCZ"),
  "Substance Use"                  = c("AUDIT", "ASI", "CANNABIS", "CPD", "DPW", "EVERSMOKE", "SMCESS")
)

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

ordered_phenotypes <- unique(unlist(phenotype_categories))

## 4. Load WLS direct and population effects ----------------------------------

direct_data <- readr::read_delim(
  direct_file,
  delim     = "\t",
  col_types = cols()
)

pop_data <- readr::read_delim(
  pop_file,
  delim     = "\t",
  col_types = cols()
)

# Merge in phenotype names
direct_data <- direct_data %>%
  dplyr::left_join(pheno_mapping, by = "Phenotype")

pop_data <- pop_data %>%
  dplyr::left_join(pheno_mapping, by = "Phenotype")

# Merge direct and population effects by phenotype
#   - direct_data: causal effects (mixed models in relatives)
#   - pop_data   : population associations (GLM/LM in unrelated sample)
merged_data <- dplyr::inner_join(
  direct_data,
  pop_data,
  by = "Phenotype",
  suffix = c("_Direct", "_Pop")
)

# Identify binary phenotypes by presence of an odds ratio for the proband effect
merged_data <- merged_data %>%
  dplyr::mutate(
    IsBinary = !is.na(OddsRatio_Proband)
  )

binary_data     <- merged_data %>% dplyr::filter(IsBinary)
non_binary_data <- merged_data %>% dplyr::filter(!IsBinary)

## 5. Panel A – Binary traits (odds ratios) ------------------------------------
#   - EffectType = "Causal Effect": OddsRatio_Proband (from GLMM in relatives)
#   - EffectType = "Population Association": OddsRatio_PopEffect_Pop (from GLM in unrelateds)

binary_data_long <- binary_data %>%
  dplyr::select(
    Phenotype,
    FullPhenoName = FullPhenoName_Direct,
    OddsRatio_Proband,
    SE_Proband,
    OddsRatio_PopEffect_Pop = OddsRatio_PopEffect_Pop,
    SE_PopEffect_Pop        = SE_PopEffect_Pop
  ) %>%
  tidyr::pivot_longer(
    cols      = c(OddsRatio_Proband, OddsRatio_PopEffect_Pop),
    names_to  = "EffectType",
    values_to = "Effect"
  ) %>%
  dplyr::mutate(
    SE = dplyr::if_else(
      EffectType == "OddsRatio_Proband",
      SE_Proband,          # SE of log-OR from GLMM
      SE_PopEffect_Pop     # SE of log-OR from GLM
    ),
    EffectType = dplyr::if_else(
      EffectType == "OddsRatio_Proband",
      "Causal Effect",
      "Population Association"
    )
  )

# 95% CIs on odds ratio scale (SEs are on log-odds scale)
binary_data_long <- binary_data_long %>%
  dplyr::mutate(
    LogEffect    = log(Effect),
    Lower_CI_log = LogEffect - 1.96 * SE,
    Upper_CI_log = LogEffect + 1.96 * SE,
    Lower_CI     = exp(Lower_CI_log),
    Upper_CI     = exp(Upper_CI_log)
  )

# Ordering by phenotype categories
binary_data_long <- binary_data_long %>%
  dplyr::mutate(
    Phenotype = factor(Phenotype, levels = ordered_phenotypes)
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

# Simple y-range and breaks including OR=1
y_min <- floor(min(binary_data_long$Lower_CI, na.rm = TRUE) * 2) / 2
y_max <- ceiling(max(binary_data_long$Upper_CI, na.rm = TRUE) * 2) / 2
if (y_min < 0.1) y_min <- 0.1
y_breaks <- pretty(c(y_min, y_max))
if (!1 %in% y_breaks) y_breaks <- sort(c(y_breaks, 1))

p_figS4A <- ggplot2::ggplot(
  binary_data_long,
  ggplot2::aes(x = FullPhenoName, y = Effect, color = EffectType, group = EffectType)
) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.6),
    size     = 3
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Lower_CI, ymax = Upper_CI),
    position = ggplot2::position_dodge(width = 0.6),
    width    = 0.3,
    linewidth = 0.6
  ) +
  ggplot2::scale_color_manual(
    values = c("Causal Effect" = "darkblue", "Population Association" = "goldenrod")
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks,
    limits = c(min(y_breaks), max(y_breaks)),
    expand = ggplot2::expansion(mult = c(0.05, 0.1))
  ) +
  ggplot2::labs(
    x     = NULL,
    y     = "Odds Ratio",
    color = NULL
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales = "free_x",
    space  = "free_x"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.background      = ggplot2::element_rect(fill = "white", color = "white"),
    panel.background     = ggplot2::element_rect(fill = "white", color = "white"),
    panel.border         = ggplot2::element_rect(colour = "grey80", fill = NA, size = 0.7),
    panel.grid.major.x   = ggplot2::element_line(color = "grey90", size = 0.3),
    panel.grid.major.y   = ggplot2::element_line(color = "grey90", size = 0.3),
    panel.grid.minor     = ggplot2::element_blank(),
    axis.text.x          = ggplot2::element_text(angle = 60, hjust = 1),
    strip.placement      = "outside",
    strip.background     = ggplot2::element_blank(),
    strip.text.x         = ggplot2::element_text(face = "bold"),
    legend.position      = "bottom",
    legend.title         = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = figS4A_file,
  plot     = p_figS4A,
  width    = 14,
  height   = 8
)

cat("Figure S4 Panel A (WLS binary causal vs population) saved to:\n", figS4A_file, "\n")

## 6. Panel B – Quantitative traits (standardized coefficients) ----------------
#   - EffectType = "Causal Effect"         : Coef_Proband (from LMM in relatives)
#   - EffectType = "Population Association": Coef_PopEffect_Pop (from LM in unrelateds)

non_binary_data_long <- non_binary_data %>%
  dplyr::select(
    Phenotype,
    FullPhenoName = FullPhenoName_Direct,
    Coef_Proband,
    SE_Proband,
    Coef_PopEffect_Pop = Coef_PopEffect_Pop,
    SE_PopEffect_Pop   = SE_PopEffect_Pop
  ) %>%
  tidyr::pivot_longer(
    cols      = c(Coef_Proband, Coef_PopEffect_Pop),
    names_to  = "EffectType",
    values_to = "Effect"
  ) %>%
  dplyr::mutate(
    SE = dplyr::if_else(
      EffectType == "Coef_Proband",
      SE_Proband,
      SE_PopEffect_Pop
    ),
    EffectType = dplyr::if_else(
      EffectType == "Coef_Proband",
      "Causal Effect",
      "Population Association"
    )
  )

# Ordering of traits and full names
non_binary_data_long <- non_binary_data_long %>%
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

# Optional: exclude specific phenotypes to match manuscript presentation
non_binary_data_long <- non_binary_data_long %>%
  dplyr::filter(
    !Phenotype %in% c("BL_HDL", "BL_TRYG", "BPpulse", "BL_CHOL", "CAD", "HARDCAD")
  )

# Category assignment
non_binary_data_long <- non_binary_data_long %>%
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

# Y-axis limits and breaks for standardized effects
y_min_nb <- floor(min(non_binary_data_long$Effect - 1.96 * non_binary_data_long$SE, na.rm = TRUE) * 2) / 2
y_max_nb <- ceiling(max(non_binary_data_long$Effect + 1.96 * non_binary_data_long$SE, na.rm = TRUE) * 2) / 2
y_breaks_nb <- pretty(c(y_min_nb, y_max_nb))

p_figS4B <- ggplot2::ggplot(
  non_binary_data_long,
  ggplot2::aes(x = FullPhenoName, y = Effect, color = EffectType, group = EffectType)
) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.6),
    size     = 3
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = Effect - 1.96 * SE,
                 ymax = Effect + 1.96 * SE),
    position = ggplot2::position_dodge(width = 0.6),
    width    = 0.3,
    linewidth = 0.6
  ) +
  ggplot2::scale_color_manual(
    values = c("Causal Effect" = "darkblue", "Population Association" = "goldenrod")
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks_nb,
    expand = ggplot2::expansion(mult = c(0.05, 0.1))
  ) +
  ggplot2::labs(
    x     = NULL,
    y     = "Standardized Effect Size",
    color = NULL
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales = "free_x",
    space  = "free_x"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    panel.background     = ggplot2::element_rect(fill = "white"),
    plot.background      = ggplot2::element_rect(fill = "white", color = "white"),
    panel.border         = ggplot2::element_rect(colour = "grey80", fill = NA, size = 0.7),
    axis.text.x          = ggplot2::element_text(angle = 60, hjust = 1),
    axis.title.y         = ggplot2::element_text(),
    legend.position      = "none",
    strip.placement      = "outside",
    strip.background     = ggplot2::element_blank(),
    strip.text.x         = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = figS4B_file,
  plot     = p_figS4B,
  width    = 14,
  height   = 8
)

cat("Figure S4 Panel B (WLS quantitative causal vs population) saved to:\n", figS4B_file, "\n")