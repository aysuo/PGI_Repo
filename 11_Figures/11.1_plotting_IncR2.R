#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------#
# Plot incremental R2 results across UKB3, HRS, and WLS.
#
# Arguments:
#   1. UKB3 results file
#   2. HRS results file
#   3. WLS results file
#   4. Full phenotype-name lookup file
#   5. Output directory
#   6. Optional personal R-library directory; use NA to omit
#
# Example:
# Rscript 11.1_plotting_IncR2_refactored.R \
#   /path/to/UKB3_EUR_r2.txt \
#   /path/to/HRS_EUR_r2.txt \
#   /path/to/WLS_EUR_r2.txt \
#   /path/to/fullPheno_name_list.txt \
#   /path/to/output \
#   NA
#----------------------------------------------------------------------------------#

########################################################
######################## Set-up ########################
########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: Rscript script.R <UKB3_r2_file> <HRS_r2_file> <WLS_r2_file> ",
    "<full_pheno_names_file> <output_dir> [library_dir|NA]"
  )
}

ukb_path <- normalizePath(args[1], mustWork = TRUE)
hrs_path <- normalizePath(args[2], mustWork = TRUE)
wls_path <- normalizePath(args[3], mustWork = TRUE)
full_names_path <- normalizePath(args[4], mustWork = TRUE)
output_dir <- normalizePath(args[5], mustWork = FALSE)

library_dir <- if (
  length(args) >= 6 &&
  nzchar(args[6]) &&
  toupper(args[6]) != "NA"
) {
  normalizePath(args[6], mustWork = FALSE)
} else {
  NULL
}

if (!is.null(library_dir)) {
  dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(library_dir, .libPaths()))
}

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "grid",
  "gridExtra",
  "cowplot"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(cowplot)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(output_dir, "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

########################################################
###################### Load data #######################
########################################################

read_result_file <- function(path, cohort_name) {
  data.table::fread(path) |>
    as_tibble() |>
    mutate(cohort = cohort_name)
}

ukb_data <- read_result_file(ukb_path, "UKB3")
hrs_data <- read_result_file(hrs_path, "HRS")
wls_data <- read_result_file(wls_path, "WLS")

combined_data <- bind_rows(ukb_data, hrs_data, wls_data)

required_result_columns <- c(
  "phenotype",
  "r2_inc_single_SBayesR",
  "r2_inc_single_SBayesR_lower",
  "r2_inc_single_SBayesR_upper"
)

# Backward compatibility with the older column naming convention.
old_to_new_names <- c(
  r2_inc_single_SBayesRlower = "r2_inc_single_SBayesR_lower",
  r2_inc_single_SBayesRupper = "r2_inc_single_SBayesR_upper"
)

for (old_name in names(old_to_new_names)) {
  new_name <- old_to_new_names[[old_name]]
  if (old_name %in% names(combined_data) && !new_name %in% names(combined_data)) {
    combined_data <- combined_data |>
      rename(!!new_name := all_of(old_name))
  }
}

missing_result_columns <- setdiff(required_result_columns, names(combined_data))

if (length(missing_result_columns) > 0) {
  stop(
    "Result files are missing required columns: ",
    paste(missing_result_columns, collapse = ", ")
  )
}

cohort_levels <- c("UKB3", "HRS", "WLS")

all_combinations <- tidyr::expand_grid(
  phenotype = unique(combined_data$phenotype),
  cohort = cohort_levels
)

combined_data <- all_combinations |>
  left_join(combined_data, by = c("phenotype", "cohort")) |>
  mutate(cohort = factor(cohort, levels = cohort_levels))

data.table::fwrite(
  combined_data,
  file.path(output_dir, "combined_data.txt"),
  sep = "\t"
)

full_pheno_names <- data.table::fread(
  full_names_path,
  header = FALSE,
  sep = ",",
  col.names = c("phenotype", "full_name")
) |>
  as_tibble() |>
  mutate(
    phenotype = gsub('^"|"$', "", phenotype),
    full_name = gsub('^"|"$', "", full_name),
    full_name = gsub(" ", "\n", full_name)
  )

duplicated_pheno <- unique(
  full_pheno_names$phenotype[duplicated(full_pheno_names$phenotype)]
)

if (length(duplicated_pheno) > 0) {
  stop(
    "Duplicate phenotype labels in full-name lookup: ",
    paste(duplicated_pheno, collapse = ", ")
  )
}

########################################################
################ Plot configuration ####################
########################################################

categories <- list(
  Anthropometric = c("BMI", "HEIGHT"),
  Blood_Biomarkers = c(
    "BL_HDL", "BL_LDL", "BL_nonHDL", "BL_CHOL",
    "BL_TRYG", "BPdia", "BPpulse", "BPsys"
  ),
  Cognition_Education = c("CP", "EA", "ALZnoproxy", "ALZ"),
  Personality_Wellbeing = c(
    "ADVENTURE", "EXTRA", "FAMSAT", "FRIENDSAT", "MORNING",
    "NEURO", "OPEN", "ACTIVITY", "RELIGATT", "RISK", "SWB"
  ),
  Health_HealthBehaviors = c(
    "ALLERGYPOLLEN", "ASTHMA", "ASTECZRHI", "BRCA", "COPD",
    "CAD", "HAYFEVER", "IBD", "MIGRAINE", "NEARSIGHTED",
    "PRCA", "SELFHEALTH", "T2D"
  ),
  Fertility_Sexual_Development = c("AFB", "MENARCHE", "AFS", "NEBwomen"),
  Psychiatric_Conditions = c(
    "ANOREX", "ADHD", "ASD", "BIPOLAR", "DEP", "INSOMNIA", "SCZ"
  ),
  Substance_Use = c(
    "AUDIT", "ASI", "CANNABIS", "CPD", "DPW", "EVERSMOKE", "SMCESS"
  )
)

category_labels <- c(
  Anthropometric = "Anthropometric",
  Blood_Biomarkers = "Blood Biomarkers",
  Cognition_Education = "Cognition & Education",
  Personality_Wellbeing = "Personality & Wellbeing",
  Health_HealthBehaviors = "Health",
  Fertility_Sexual_Development = "Fertility & Sexual Development",
  Psychiatric_Conditions = "Psychiatric Conditions",
  Substance_Use = "Substance Use"
)

colour_scale <- c(
  UKB3 = "#33a02c",
  HRS = "#1f78b4",
  WLS = "#ff7f00"
)

ordered_phenotypes <- unlist(categories, use.names = FALSE)

missing_full_names <- setdiff(
  intersect(ordered_phenotypes, combined_data$phenotype),
  full_pheno_names$phenotype
)

if (length(missing_full_names) > 0) {
  stop(
    "Missing full names for phenotypes present in the results: ",
    paste(missing_full_names, collapse = ", ")
  )
}

plot_data <- combined_data |>
  mutate(
    pheno_group = case_when(
      phenotype %in% categories$Anthropometric ~ "Anthropometric",
      phenotype %in% categories$Blood_Biomarkers ~ "Blood_Biomarkers",
      phenotype %in% categories$Cognition_Education ~ "Cognition_Education",
      phenotype %in% categories$Personality_Wellbeing ~ "Personality_Wellbeing",
      phenotype %in% categories$Health_HealthBehaviors ~ "Health_HealthBehaviors",
      phenotype %in% categories$Fertility_Sexual_Development ~ "Fertility_Sexual_Development",
      phenotype %in% categories$Psychiatric_Conditions ~ "Psychiatric_Conditions",
      phenotype %in% categories$Substance_Use ~ "Substance_Use",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(pheno_group)) |>
  left_join(full_pheno_names, by = "phenotype") |>
  mutate(
    pheno_group_fullname = unname(category_labels[pheno_group]),
    phenotype = factor(phenotype, levels = ordered_phenotypes),
    full_name = factor(
      full_name,
      levels = full_pheno_names |>
        filter(phenotype %in% ordered_phenotypes) |>
        arrange(match(phenotype, ordered_phenotypes)) |>
        pull(full_name)
    )
  )

########################################################
############ Individual category barplots ##############
########################################################

plot_category <- function(category_name, phenotypes) {
  category_data <- combined_data |>
    filter(phenotype %in% phenotypes) |>
    complete(
      phenotype,
      cohort = factor(cohort_levels, levels = cohort_levels)
    ) |>
    mutate(
      phenotype = factor(phenotype, levels = phenotypes)
    )

  finite_upper <- category_data$r2_inc_single_SBayesR_upper[
    is.finite(category_data$r2_inc_single_SBayesR_upper)
  ]

  if (length(finite_upper) == 0) {
    warning("No finite values for category: ", category_name)
    return(NULL)
  }

  ymax <- max(1, ceiling(max(finite_upper)))
  break_size <- max(1, ceiling(ymax / 5))

  plot <- ggplot(
    category_data,
    aes(
      x = phenotype,
      y = r2_inc_single_SBayesR,
      fill = cohort
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.75,
      colour = "black",
      na.rm = TRUE
    ) +
    geom_errorbar(
      aes(
        ymin = r2_inc_single_SBayesR_lower,
        ymax = r2_inc_single_SBayesR_upper
      ),
      position = position_dodge(width = 0.8),
      width = 0.15,
      na.rm = TRUE
    ) +
    scale_fill_manual(values = colour_scale, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(0, ymax, by = break_size),
      limits = c(0, ymax),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = category_labels[[category_name]],
      y = expression(Incremental ~ italic(R^2))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(category_name, "_barplot.png")
    ),
    plot = plot,
    width = 10,
    height = 7,
    dpi = 300
  )

  plot
}

category_plots <- purrr::imap(
  categories,
  ~ plot_category(.y, .x)
)

category_plots <- category_plots[
  !vapply(category_plots, is.null, logical(1))
]

if (length(category_plots) > 0) {
  category_legend_plot <- ggplot(
    combined_data,
    aes(
      x = phenotype,
      y = r2_inc_single_SBayesR,
      fill = cohort
    )
  ) +
    geom_col() +
    scale_fill_manual(values = colour_scale, drop = FALSE) +
    theme_void() +
    theme(legend.position = "bottom")

  category_legend <- cowplot::get_legend(category_legend_plot)

  category_grobs <- lapply(category_plots, ggplotGrob)
  composite_categories <- gridExtra::arrangeGrob(
    grobs = c(category_grobs, list(category_legend)),
    ncol = 3
  )

  ggsave(
    filename = file.path(figure_dir, "Composite_Figure.png"),
    plot = composite_categories,
    width = 20,
    height = 15,
    dpi = 300
  )
}

########################################################
################ Grouped manuscript plots ##############
########################################################

base_group_theme <- theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 32, colour = "grey4"),
    axis.text.x = element_text(
      size = 28,
      colour = "grey12",
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(size = 20, colour = "grey12"),
    strip.placement = "outside",
    strip.text.x = element_text(size = 32),
    strip.background = element_rect(fill = NA, colour = "white"),
    panel.spacing = unit(0, "lines"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(
      colour = "#bcacac",
      fill = NA,
      linewidth = 0.5
    )
  )

make_grouped_plot <- function(
  groups,
  y_breaks,
  y_limit = NULL,
  show_legend = FALSE,
  exclusions = NULL
) {
  current_data <- plot_data |>
    filter(
      pheno_group %in% groups,
      is.finite(r2_inc_single_SBayesR)
    )

  if (!is.null(exclusions)) {
    current_data <- exclusions(current_data)
  }

  plot <- current_data |>
    ggplot(
      aes(
        fill = cohort,
        x = full_name,
        y = r2_inc_single_SBayesR
      )
    ) +
    facet_grid(
      . ~ pheno_group_fullname,
      space = "free",
      scales = "free",
      switch = "x"
    ) +
    geom_col(
      colour = "black",
      width = 0.6,
      position = position_dodge2(
        preserve = "single",
        padding = 0
      )
    ) +
    geom_errorbar(
      aes(
        ymin = r2_inc_single_SBayesR_lower,
        ymax = r2_inc_single_SBayesR_upper
      ),
      width = 0.6,
      position = position_dodge2(
        width = 0,
        preserve = "single",
        padding = 0.5
      )
    ) +
    base_group_theme +
    scale_fill_manual(values = colour_scale, drop = FALSE) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = y_limit,
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.07))
    ) +
    labs(
      x = NULL,
      y = expression(Incremental ~ italic(R^2))
    ) +
    coord_cartesian(clip = "off")

  if (show_legend) {
    plot <- plot +
      theme(
        legend.title = element_blank(),
        legend.position = c(0.96, 0.8),
        legend.background = element_rect(
          colour = "black",
          fill = "white"
        )
      )
  } else {
    plot <- plot + theme(legend.position = "none")
  }

  plot
}

plot123 <- make_grouped_plot(
  groups = c(
    "Anthropometric",
    "Blood_Biomarkers",
    "Cognition_Education"
  ),
  y_breaks = seq(0, 45, 5),
  show_legend = TRUE,
  exclusions = function(df) {
    df |>
      filter(
        !(
          cohort %in% c("WLS", "HRS") &
            phenotype %in% c(
              "BL_HDL",
              "BL_LDL",
              "BL_TRYG",
              "BL_CHOL"
            )
        ),
        !(cohort == "WLS" & phenotype == "BPpulse")
      )
  }
)

plot4 <- make_grouped_plot(
  groups = c(
    "Health_HealthBehaviors",
    "Fertility_Sexual_Development"
  ),
  y_breaks = seq(0, 13, 2),
  y_limit = c(0, 13)
)

plot5 <- make_grouped_plot(
  groups = "Personality_Wellbeing",
  y_breaks = seq(0, 6, 2),
  y_limit = c(0, 6)
)

plot6 <- make_grouped_plot(
  groups = c("Psychiatric_Conditions", "Substance_Use"),
  y_breaks = seq(0, 7, 2),
  y_limit = c(0, 7)
)

grouped_plots <- list(
  plot123 = plot123,
  plot4 = plot4,
  plot5 = plot5,
  plot6 = plot6
)

for (plot_name in names(grouped_plots)) {
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(plot_name, ".png")
    ),
    plot = grouped_plots[[plot_name]],
    width = 10,
    height = 6,
    dpi = 300
  )
}

combined_group_plot <- gridExtra::arrangeGrob(
  grobs = lapply(grouped_plots, ggplotGrob),
  ncol = 1
)

ggsave(
  filename = file.path(
    figure_dir,
    "combined_plot_6phenoGroups_fullPhenoNames.png"
  ),
  plot = combined_group_plot,
  width = 36,
  height = 40,
  units = "in",
  dpi = 400,
  limitsize = FALSE
)
