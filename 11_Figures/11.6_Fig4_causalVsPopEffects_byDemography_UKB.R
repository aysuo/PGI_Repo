#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 10.6.1_Fig4_casualVsPopEffects_byDemography_UKB.R
#
# Purpose:
#   Generate Figure 4: standardized regression coefficients of PGIs for
#   diastolic blood pressure (BPdia) and BMI from models that do and do not
#   control for parental PGIs in the first-degree-relatives subsample of UKB1.
#
# Description:
#   - Inputs are phenotype- and interaction-specific coefficient files created
#     by the upstream mixed model scripts:
#       • BPdia: sex × PGI interaction
#       • BMI  : age-group × PGI interaction (below vs above median age)
#
#   - For each phenotype, the script:
#       1) Reconstructs the PGI effect in the baseline group
#          (Male / Younger) and the comparison group (Female / Older),
#          separately for:
#             • Population Association model (no parental PGIs)
#             • Causal Effect model (adjusted for parental PGIs)
#       2) Computes 95% confidence intervals for each effect.
#
#   - Outputs a single two-panel figure:
#       • Panel A: PGI effect on diastolic blood pressure by sex.
#       • Panel B: PGI effect on BMI by age group (below vs above median).
#
#   Shapes:
#       • Circle  (●): models that do NOT control for parental PGIs
#                      (“Population Association”).
#       • Triangle (▲): models that DO control for parental PGIs
#                       (“Causal Effect”).
#
#   Error bars represent 95% confidence intervals.
#
# Environment variables (set in paths_PGIrepo_withinFam):
#   - PGI_RepoV2               : base project directory
#   - UKB1_DEMOG_RESULTS_DIR   : directory containing
#         PGI_by_Sex_interaction_results_for_plotting_BPdia.csv
#         PGI_by_Sex_interaction_results_for_plotting_BMI.csv
#         (default: <PGI_RepoV2>/processed/sumStats/within_family_related)
#   - UKB1_FIG_DIR             : output directory for the figure
#         (default: <PGI_RepoV2>/processed/7_Figures/direct_indirect_PGIs)
# ------------------------------------------------------------------------------

# ── 1. Setup paths ────────────────────────────────────────────────────────────

PGI_RepoV2 <- Sys.getenv("PGI_RepoV2")
if (PGI_RepoV2 == "") {
  stop("PGI_RepoV2 is not defined. Please set it in your environment before running.")
}

demog_dir_env <- Sys.getenv("UKB1_DEMOG_RESULTS_DIR")
if (demog_dir_env == "") {
  demog_dir <- file.path(
    PGI_RepoV2,
    "processed", "sumStats", "within_family_related"
  )
} else {
  demog_dir <- demog_dir_env
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

bpdia_file <- file.path(
  demog_dir,
  "PGI_by_Sex_interaction_results_for_plotting_BPdia.csv"
)
bmi_file <- file.path(
  demog_dir,
  "PGI_by_Sex_interaction_results_for_plotting_BMI.csv"
)

output_file_full <- file.path(
  fig_dir,
  "Full_Effect_Comparison_Plot_BPdia_BMI.png"
)

# ── 2. Libraries ──────────────────────────────────────────────────────────────

personal_lib_path <- file.path(PGI_RepoV2, "doc", "R", "library")
.libPaths(c(personal_lib_path, .libPaths()))

pkgs <- c("dplyr", "readr", "tidyr", "ggplot2", "purrr", "stringr")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = personal_lib_path)) {
    install.packages(pkg, lib = personal_lib_path, repos = "https://cran.r-project.org/")
    library(pkg, character.only = TRUE, lib.loc = personal_lib_path)
  }
}

# ── 3. Phenotype-specific metadata ────────────────────────────────────────────

phenotype_info <- list(
  BPdia = list(
    file_path      = bpdia_file,
    interaction_term = "proband:SexFemale",
    panel_title    = "A. Diastolic Blood Pressure",
    baseline_label = "Male",
    comparison_label = "Female"
  ),
  BMI = list(
    file_path      = bmi_file,
    interaction_term = "proband:AgeGroupOlder",
    panel_title    = "B. Body Mass Index",
    baseline_label = "Younger",
    comparison_label = "Older"
  )
)

# ── 4. Helper: prepare data for plotting for one phenotype ────────────────────
#   Input:  one phenotype’s coefficient CSV produced upstream.
#   Output: long-format tibble with:
#           phenotype, panel_title, model, group, estimate, lower, upper

prepare_data_for_plotting <- function(pheno_name, info) {

  if (!file.exists(info$file_path)) {
    warning("File not found for ", pheno_name, ": ", info$file_path)
    return(NULL)
  }

  raw_data <- readr::read_csv(info$file_path, show_col_types = FALSE) %>%
    dplyr::filter(phenotype == pheno_name)

  if (nrow(raw_data) == 0) {
    warning("No data rows for phenotype '", pheno_name,
            "' in file: ", info$file_path)
    return(NULL)
  }

  # We need the proband and interaction terms, by model (pop/direct)
  wide_betas <- raw_data %>%
    dplyr::filter(term %in% c("proband", info$interaction_term)) %>%
    tidyr::pivot_wider(
      id_cols    = c(phenotype, model),
      names_from = term,
      values_from = c(estimate, std.error, covariance)
    )

  required_cols <- c("estimate_proband",
                     paste0("estimate_", info$interaction_term))
  if (!all(required_cols %in% names(wide_betas))) {
    warning("Missing required beta terms for phenotype '", pheno_name, "'.")
    return(NULL)
  }

  # Extract incremental R² is not needed for this figure, so ignore it.
  processed <- wide_betas %>%
    dplyr::rename(
      beta_baseline    = estimate_proband,
      se_baseline      = std.error_proband,
      beta_interaction = !!rlang::sym(paste0("estimate_", info$interaction_term)),
      se_interaction   = !!rlang::sym(paste0("std.error_", info$interaction_term)),
      cov_bi           = !!rlang::sym(paste0("covariance_", info$interaction_term))
    ) %>%
    dplyr::mutate(
      cov_bi   = dplyr::if_else(is.na(cov_bi), 0, cov_bi),
      # Comparison (Female / Older) = baseline + interaction
      beta_comparison = beta_baseline + beta_interaction,

      var_baseline    = se_baseline^2,
      var_interaction = se_interaction^2,
      var_comparison  = var_baseline + var_interaction + 2 * cov_bi,

      baseline_low  = beta_baseline  - 1.96 * se_baseline,
      baseline_high = beta_baseline  + 1.96 * se_baseline,
      comparison_low  = beta_comparison - 1.96 * sqrt(pmax(var_comparison, 0)),
      comparison_high = beta_comparison + 1.96 * sqrt(pmax(var_comparison, 0))
    )

  # Long format: rows = (model × group)
  out <- processed %>%
    dplyr::select(
      phenotype, model,
      beta_baseline, baseline_low, baseline_high,
      beta_comparison, comparison_low, comparison_high
    ) %>%
    tidyr::pivot_longer(
      cols = c(beta_baseline, beta_comparison),
      names_to = "group_type",
      values_to = "estimate"
    ) %>%
    dplyr::mutate(
      lower = dplyr::case_when(
        group_type == "beta_baseline"   ~ baseline_low,
        group_type == "beta_comparison" ~ comparison_low
      ),
      upper = dplyr::case_when(
        group_type == "beta_baseline"   ~ baseline_high,
        group_type == "beta_comparison" ~ comparison_high
      ),
      group = dplyr::case_when(
        group_type == "beta_baseline"   ~ info$baseline_label,
        group_type == "beta_comparison" ~ info$comparison_label
      ),
      panel_title = info$panel_title
    ) %>%
    dplyr::select(
      phenotype, panel_title, model, group, estimate, lower, upper
    )

  return(out)
}

# ── 5. Build combined plotting dataset for BPdia and BMI ──────────────────────

processed_data <- purrr::imap_dfr(
  phenotype_info,
  ~ prepare_data_for_plotting(pheno_name = .y, info = .x)
)

if (is.null(processed_data) || nrow(processed_data) == 0) {
  stop("No data available for plotting. Check that the input CSV files exist and are non-empty.")
}

# Recode model labels and set factor levels
plot_data <- processed_data %>%
  dplyr::mutate(
    model = dplyr::recode(
      model,
      "pop"    = "Population Association",
      "direct" = "Causal Effect"
    ),
    model = factor(model, levels = c("Population Association", "Causal Effect")),
    group = factor(group, levels = c("Male", "Female", "Younger", "Older")),
    panel_title = factor(
      panel_title,
      levels = c("A. Diastolic Blood Pressure", "B. Body Mass Index")
    )
  )

# ── 6. Generate Figure 4: full effect comparison plot ─────────────────────────

# Mapping shapes: circle = no parental PGIs; triangle = with parental PGIs
shape_values <- c(
  "Population Association" = 16,  # circle
  "Causal Effect"          = 17   # triangle
)

full_comparison_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = group,
    y = estimate,
    shape = model
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    color      = "grey40"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    width     = 0.2,
    linewidth = 0.8,
    position  = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::geom_point(
    size     = 3,
    position = ggplot2::position_dodge(width = 0.5)
  ) +
  ggplot2::facet_wrap(
    ~ panel_title,
    ncol   = 1,
    scales = "free_x"
  ) +
  ggplot2::scale_shape_manual(
    values = shape_values,
    name   = "Model"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Standardized PGI Effect (95% CI)"
  ) +
  ggplot2::theme_bw(base_size = 14) +
  ggplot2::theme(
    legend.position   = "bottom",
    legend.title      = ggplot2::element_blank(),
    strip.background  = ggplot2::element_rect(fill = "grey90", color = "grey50"),
    strip.text        = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
    axis.text.x       = ggplot2::element_text(size = 11),
    axis.text.y       = ggplot2::element_text(size = 11),
    plot.title        = ggplot2::element_blank(),
    plot.subtitle     = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = output_file_full,
  plot     = full_comparison_plot,
  width    = 8,
  height   = 9,
  dpi      = 300
)

cat("\n✔ Full effect comparison plot for BPdia and BMI saved to:\n  ",
    output_file_full, "\n")