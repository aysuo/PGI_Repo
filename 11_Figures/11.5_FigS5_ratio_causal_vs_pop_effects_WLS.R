#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.5.4_FigS5_ratio_causal_vs_pop_effects_WLS.R
# Purpose: Recreate Figure S5 showing the ratio of causal effects of PGIs to
#          their population associations in the WLS cohort, for both binary
#          and quantitative phenotypes.
#
# Description:
#   - Causal effects:
#       • Binary traits: log-odds effects from generalized linear mixed models
#         with a logistic link and a random intercept for family ID, estimated
#         in first-degree relatives (sibling genotyped sample).
#       • Quantitative traits: standardized coefficients from linear mixed
#         models with a random intercept for family ID in the relative sample.
#
#   - Population associations:
#       • Binary traits: log-odds from logistic regression in the sample of
#         unrelated individuals.
#       • Quantitative traits: standardized coefficients from linear regression
#         in the unrelated sample.
#
#   - Ratios:
#       • Binary traits: ratio of log-odds effects
#             Ratio = log(OR_causal) / log(OR_population)
#       • Quantitative traits: ratio of standardized coefficients
#             Ratio = beta_causal / beta_population
#       • Standard errors of ratios are computed using the delta method.
#
#   Panels:
#     • Panel A: Binary traits, ratio on the log-effect scale (y-axis shows
#                Ratio = logOR_causal / logOR_pop).
#     • Panel B: Quantitative traits, ratio of standardized linear effects.
#
# Inputs (created by previous WLS scripts):
#   - <WLS_OUTPUT_DIR>/binaryNonBinary_directPGIeffect_coefficients.txt
#       (causal effects in relatives; from 10.4.1_causalEffect_WLS_BinaryTraits_lme4.R)
#   - <WLS_OUTPUT_DIR>/binaryNonBinary_popPGIeffect_coefficients.txt
#       (population associations in unrelateds; from 10.4.2_populationAssociation_WLS_allTraits.R)
#   - <PHENO_MAPPING_FILE> (tab-delimited: Phenotype, FullPhenoName)
#
# Outputs:
#   - <WLS_FIG_DIR>/FigS5A_WLS_binary_ratio_causal_vs_population.png
#   - <WLS_FIG_DIR>/FigS5B_WLS_quantitative_ratio_causal_vs_population.png
#   - <WLS_SUMSTAT_DIR>/summaryStat_binaryPheno_WLS.csv
#   - <WLS_SUMSTAT_DIR>/summaryStat_nonbinaryPheno_WLS.csv
#
# Environment variables (e.g. set via paths_PGIrepo_withinFam):
#   - PGI_RepoV2         : base project directory
#   - WLS_OUTPUT_DIR     : directory with WLS direct & population effect files
#   - WLS_FIG_DIR        : directory for WLS figures (optional; defaults to
#                          <PGI_RepoV2>/figures if unset)
#   - WLS_SUMSTAT_DIR    : directory for WLS summary statistics (optional;
#                          defaults to <PGI_RepoV2>/processed/sumStats/WLS)
#   - PHENO_MAPPING_FILE : file with short ↔ full phenotype names
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
wls_sumstat_dir_env <- Sys.getenv("WLS_SUMSTAT_DIR")
pheno_mapping_file <- Sys.getenv("PHENO_MAPPING_FILE")

if (wls_output_dir == "") {
  stop("WLS_OUTPUT_DIR is not defined. Set it in your environment (paths_PGIrepo_withinFam).")
}
if (pheno_mapping_file == "") {
  stop("PHENO_MAPPING_FILE is not defined. Set it in your environment (paths_PGIrepo_withinFam).")
}

# Figure directory
if (wls_fig_dir_env == "") {
  wls_fig_dir <- file.path(PGI_RepoV2, "figures")
} else {
  wls_fig_dir <- wls_fig_dir_env
}
if (!dir.exists(wls_fig_dir)) {
  dir.create(wls_fig_dir, recursive = TRUE, showWarnings = FALSE)
}

figS5A_file <- file.path(wls_fig_dir, "FigS5A_WLS_binary_ratio_causal_vs_population.png")
figS5B_file <- file.path(wls_fig_dir, "FigS5B_WLS_quantitative_ratio_causal_vs_population.png")

# Summary-stat directory
if (wls_sumstat_dir_env == "") {
  wls_sumstat_dir <- file.path(PGI_RepoV2, "processed", "sumStats", "WLS")
} else {
  wls_sumstat_dir <- wls_sumstat_dir_env
}
if (!dir.exists(wls_sumstat_dir)) {
  dir.create(wls_sumstat_dir, recursive = TRUE, showWarnings = FALSE)
}

binary_summary_file     <- file.path(wls_sumstat_dir, "summaryStat_binaryPheno_WLS.csv")
nonbinary_summary_file  <- file.path(wls_sumstat_dir, "summaryStat_nonbinaryPheno_WLS.csv")

# Direct & population effect files
direct_file <- file.path(
  wls_output_dir,
  "binaryNonBinary_directPGIeffect_coefficients.txt"
)
pop_file <- file.path(
  wls_output_dir,
  "binaryNonBinary_popPGIeffect_coefficients.txt"
)

## 2. Libraries ----------------------------------------------------------------

personal_lib_path <- file.path(PGI_RepoV2, "doc", "R", "library")
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

pheno_mapping <- readr::read_delim(
  file      = pheno_mapping_file,
  delim     = "\t",
  col_names = c("Phenotype", "FullPhenoName"),
  col_types = cols()
)

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

## 4. Load and merge WLS direct & population effects --------------------------

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

direct_data <- direct_data %>%
  dplyr::left_join(pheno_mapping, by = "Phenotype")

pop_data <- pop_data %>%
  dplyr::left_join(pheno_mapping, by = "Phenotype")

merged_data <- dplyr::inner_join(
  direct_data,
  pop_data,
  by = "Phenotype",
  suffix = c("_Direct", "_Pop")
)

# Prefer the "direct" full name column
if (!"FullPhenoName_Direct" %in% names(merged_data) &&
    "FullPhenoName_Pop" %in% names(merged_data)) {
  merged_data$FullPhenoName_Direct <- merged_data$FullPhenoName_Pop
}

## 5. Compute ratios and SEs (delta method) ------------------------------------

# Identify binary phenotypes via presence of an odds ratio on proband effect
merged_data <- merged_data %>%
  dplyr::mutate(
    IsBinary = !is.na(OddsRatio_Proband)
  )

# Ratio definitions:
#   Binary:     Ratio = log(OR_direct) / log(OR_pop)
#   Non-binary: Ratio = beta_direct / beta_pop
merged_data <- merged_data %>%
  dplyr::mutate(
    # Transform to log-odds for binary; leave as coefficients otherwise
    LogEff_Proband = dplyr::if_else(
      IsBinary,
      log(OddsRatio_Proband),
      Coef_Proband
    ),
    LogEff_Pop = dplyr::if_else(
      IsBinary,
      log(OddsRatio_PopEffect_Pop),
      Coef_PopEffect_Pop
    ),
    # Ratio of effects
    Ratio = dplyr::if_else(
      IsBinary,
      LogEff_Proband / LogEff_Pop,
      Coef_Proband / Coef_PopEffect_Pop
    ),
    # SEs on the appropriate scale:
    SE_LogEff_Proband = SE_Proband,
    SE_LogEff_Pop     = SE_PopEffect_Pop
  ) %>%
  dplyr::mutate(
    SE_Ratio = dplyr::if_else(
      IsBinary,
      # Binary traits – delta method on log-odds ratio
      abs(Ratio) * sqrt(
        (SE_LogEff_Proband / LogEff_Proband)^2 +
          (SE_LogEff_Pop     / LogEff_Pop)^2
      ),
      # Quantitative traits – delta method for a/b
      sqrt(
        (SE_Proband^2) / (Coef_PopEffect_Pop^2) +
          (Coef_Proband^2 * SE_PopEffect_Pop^2) / (Coef_PopEffect_Pop^4)
      )
    )
  ) %>%
  dplyr::mutate(
    LowerCI_Ratio = Ratio - 1.96 * SE_Ratio,
    UpperCI_Ratio = Ratio + 1.96 * SE_Ratio
  )

binary_data     <- merged_data %>% dplyr::filter(IsBinary)
non_binary_data <- merged_data %>% dplyr::filter(!IsBinary)

## 6. Panel A – Binary traits: Ratio plots -------------------------------------

binary_data <- binary_data %>%
  dplyr::mutate(
    Phenotype = factor(Phenotype, levels = ordered_phenotypes)
  ) %>%
  dplyr::arrange(Phenotype) %>%
  dplyr::mutate(
    FullPhenoName = FullPhenoName_Direct,
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

ymin_bin <- floor(min(binary_data$LowerCI_Ratio, na.rm = TRUE))
ymax_bin <- ceiling(max(binary_data$UpperCI_Ratio, na.rm = TRUE))
y_breaks_bin <- seq(ymin_bin, ymax_bin, by = 1)
y_minor_bin  <- seq(ymin_bin, ymax_bin, by = 0.5)

p_figS5A <- ggplot2::ggplot(
  binary_data,
  ggplot2::aes(x = FullPhenoName, y = Ratio)
) +
  ggplot2::geom_point(
    color    = "black",
    size     = 3,
    position = ggplot2::position_dodge(width = 0.6)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = LowerCI_Ratio, ymax = UpperCI_Ratio),
    color     = "black",
    width     = 0.3,
    linewidth = 0.6,
    position  = ggplot2::position_dodge(width = 0.6)
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype   = "dashed",
    color      = "red"
  ) +
  ggplot2::scale_y_continuous(
    breaks       = y_breaks_bin,
    minor_breaks = y_minor_bin,
    expand       = ggplot2::expansion(mult = c(0.05, 0.1))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Causal Effect : Population Association Ratio"
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales = "free_x",
    space  = "free_x"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white"),
    plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
    panel.grid.major = ggplot2::element_line(color = "grey85", size = 0.5),
    panel.grid.minor = ggplot2::element_line(color = "grey92", size = 0.3),
    axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1),
    axis.title.y     = ggplot2::element_text(),
    legend.position  = "none",
    strip.placement  = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x     = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = figS5A_file,
  plot     = p_figS5A,
  width    = 14,
  height   = 8
)

cat("Figure S5 Panel A (WLS binary ratio causal vs population) saved to:\n", figS5A_file, "\n")

## 7. Panel B – Quantitative traits: Ratio plots -------------------------------

non_binary_data <- non_binary_data %>%
  dplyr::mutate(
    Phenotype = factor(Phenotype, levels = ordered_phenotypes)
  ) %>%
  dplyr::arrange(Phenotype) %>%
  dplyr::mutate(
    FullPhenoName = FullPhenoName_Direct,
    FullPhenoName = factor(
      FullPhenoName,
      levels = unique(FullPhenoName[order(match(Phenotype, ordered_phenotypes))])
    )
  )

# Optionally exclude the same traits as in S4 quantitative panel
non_binary_data <- non_binary_data %>%
  dplyr::filter(
    !Phenotype %in% c("BL_HDL", "BL_TRYG", "BPpulse", "BL_CHOL")
  ) %>%
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

ymin_nb <- floor(min(non_binary_data$LowerCI_Ratio, na.rm = TRUE) * 2) / 2
ymax_nb <- ceiling(max(non_binary_data$UpperCI_Ratio, na.rm = TRUE) * 2) / 2
y_breaks_nb <- seq(ymin_nb, ymax_nb, by = 0.5)

p_figS5B <- ggplot2::ggplot(
  non_binary_data %>% dplyr::filter(FullPhenoName != "Coronary Artery Disease"),
  ggplot2::aes(x = FullPhenoName, y = Ratio)
) +
  ggplot2::geom_point(
    color    = "black",
    size     = 3,
    position = ggplot2::position_dodge(width = 0.6)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = LowerCI_Ratio, ymax = UpperCI_Ratio),
    color     = "black",
    width     = 0.3,
    linewidth = 0.6,
    position  = ggplot2::position_dodge(width = 0.6)
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype   = "dashed",
    color      = "red"
  ) +
  ggplot2::scale_y_continuous(
    breaks = y_breaks_nb,
    expand = ggplot2::expansion(mult = c(0.05, 0.1))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Causal Effect : Population Association Ratio"
  ) +
  ggplot2::facet_grid(
    ~ Category,
    scales = "free_x",
    space  = "free_x"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white"),
    plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
    axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1),
    axis.title.y     = ggplot2::element_text(),
    legend.position  = "none",
    strip.placement  = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x     = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = figS5B_file,
  plot     = p_figS5B,
  width    = 14,
  height   = 8
)

cat("Figure S5 Panel B (WLS quantitative ratio causal vs population) saved to:\n", figS5B_file, "\n")

## 8. Summary statistics & BH-adjusted p-values -------------------------------

# Binary summary stats
binary_summary_stats <- binary_data %>%
  dplyr::mutate(
    PopEffect_psi      = OddsRatio_PopEffect_Pop,
    psi_SE             = SE_PopEffect_Pop,
    DirectEffect_delta = OddsRatio_Proband,
    delta_SE           = SE_Proband,
    DirectPop_Ratio    = Ratio,
    SE_ratio           = SE_Ratio,
    # Z-score under H0: ratio = 1
    z_score = (DirectPop_Ratio - 1) / SE_ratio,
    p_value = 2 * (1 - pnorm(abs(z_score)))
  ) %>%
  dplyr::select(
    Phenotype,
    PopEffect_psi,
    psi_SE,
    DirectEffect_delta,
    delta_SE,
    DirectPop_Ratio,
    SE_ratio,
    z_score,
    p_value
  )

# Non-binary summary stats
non_binary_summary_stats <- non_binary_data %>%
  dplyr::mutate(
    PopEffect_psi      = Coef_PopEffect_Pop,
    psi_SE             = SE_PopEffect_Pop,
    DirectEffect_delta = Coef_Proband,
    delta_SE           = SE_Proband,
    DirectPop_Ratio    = Ratio,
    SE_ratio           = SE_Ratio,
    z_score            = (DirectPop_Ratio - 1) / SE_ratio,
    p_value            = 2 * (1 - pnorm(abs(z_score)))
  ) %>%
  dplyr::select(
    Phenotype,
    PopEffect_psi,
    psi_SE,
    DirectEffect_delta,
    delta_SE,
    DirectPop_Ratio,
    SE_ratio,
    z_score,
    p_value
  )

# Combine and adjust p-values (BH) across all phenotypes
combined_summary_stats <- dplyr::bind_rows(
  binary_summary_stats,
  non_binary_summary_stats
) %>%
  dplyr::mutate(
    adjusted_p_value = p.adjust(p_value, method = "BH")
  )

# Split back into binary / non-binary
binary_summary_stats <- combined_summary_stats %>%
  dplyr::filter(Phenotype %in% binary_data$Phenotype)

non_binary_summary_stats <- combined_summary_stats %>%
  dplyr::filter(Phenotype %in% non_binary_data$Phenotype)

# Write summary statistic files
readr::write_csv(binary_summary_stats,    binary_summary_file)
readr::write_csv(non_binary_summary_stats, nonbinary_summary_file)

cat("Binary summary stats written to:\n",    binary_summary_file,    "\n")
cat("Non-binary summary stats written to:\n", nonbinary_summary_file, "\n")