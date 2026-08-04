#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------#
# Predicts HRS/WLS/UKB phenotypes from one or more PGI score types.
# Phenotype files must be named <PHENOTYPE>_<ANCESTRY>.pheno, e.g. BMI_AFR.pheno.
# Score files are matched as PGI_<COHORT>_<ANCESTRY>_<PHENOTYPE>-<SCORETYPE>*.
#
# Arguments:
#   1. cohort
#   2. ancestry
#        Ancestry label used in score filenames, e.g. AFR, EUR, EAS.
#   3. pheno_dir
#   4. pgi_dir
#   5. crosswalk
#        Required for HRS; use NA for other cohorts.
#   6. score_types
#        Comma-separated, e.g.:
#        single_SBayesR,single_LDpred,multi_LDpred
#   7. comparisons
#        NA for no comparisons, or semicolon-separated score-type pairs.
#        Within each pair, separate score types with a colon, e.g.:
#        single_SBayesR:single_LDpred;multi_LDpred:single_LDpred
#        Each comparison is calculated as first score type minus second score type.
#   8. output_dir
#        Optional. Defaults to derived_data/10_Prediction/output
#   9. iterations
#        Optional. Defaults to 1000.
#  10. pc_path
#        Optional. Use NA to retain the cohort-specific default.
#
# Example:
# Rscript 10.2.1_predict_phenotypes_refactored_with_ancestry.R \
#   HRS \
#   AFR \
#   /path/to/phenotypes \
#   /path/to/scores \
#   /path/to/HRS_crosswalk.csv \
#   single_SBayesR,single_LDpred \
#   single_SBayesR:single_LDpred \
#   /path/to/output \
#   1000 \
#   /path/to/HRS_PCs.eigenvec
#----------------------------------------------------------------------------------#

########################################################
######################## Set-up ########################
########################################################

required_packages <- c("data.table", "dplyr", "tidyr", "stringr", "purrr")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop(
    "Usage: Rscript script.R <cohort> <ancestry> <pheno_dir> <pgi_dir> ",
    "<crosswalk|NA> <score_types> <comparisons|NA> ",
    "[output_dir] [iterations] [pc_path|NA]"
  )
}

cohort <- args[1]
ancestry <- toupper(str_trim(args[2]))

if (!nzchar(ancestry)) {
  stop("ancestry must be a non-empty label, e.g. AFR or EUR.")
}

pheno_dir <- normalizePath(args[3], mustWork = TRUE)
pgi_dir <- normalizePath(args[4], mustWork = TRUE)

crosswalk <- if (toupper(args[5]) == "NA") {
  NA_character_
} else {
  normalizePath(args[5], mustWork = TRUE)
}

score_types <- str_split(args[6], ",", simplify = TRUE) |>
  as.character() |>
  str_trim()
score_types <- score_types[nzchar(score_types)]

if (length(score_types) == 0) {
  stop("At least one score type must be supplied.")
}

comparison_arg <- args[7]

output_dir <- if (length(args) >= 8 && nzchar(args[8])) {
  normalizePath(args[8], mustWork = FALSE)
} else {
  file.path("derived_data", "10_Prediction", "output")
}

iterations <- if (length(args) >= 9 && nzchar(args[9])) {
  as.integer(args[9])
} else {
  1000L
}

if (is.na(iterations) || iterations < 1) {
  stop("iterations must be a positive integer.")
}

pc_path_arg <- if (
  length(args) >= 10 &&
  nzchar(args[10]) &&
  toupper(args[10]) != "NA"
) {
  normalizePath(args[10], mustWork = TRUE)
} else {
  NA_character_
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
bootstrap_output_dir <- file.path(output_dir, "bootstraps")
dir.create(bootstrap_output_dir, recursive = TRUE, showWarnings = FALSE)

########################################################
############ Parse requested comparisons ###############
########################################################

parse_comparisons <- function(comparison_arg, valid_score_types) {
  if (
    length(comparison_arg) == 0 ||
    is.na(comparison_arg) ||
    !nzchar(comparison_arg) ||
    toupper(comparison_arg) == "NA"
  ) {
    return(tibble(
      score_type_1 = character(),
      score_type_2 = character(),
      comparison_name = character()
    ))
  }

  pair_strings <- str_split(comparison_arg, ";", simplify = FALSE)[[1]]
  pair_strings <- str_trim(pair_strings)
  pair_strings <- pair_strings[nzchar(pair_strings)]

  parsed <- map_dfr(pair_strings, function(pair_string) {
    pair <- str_split(pair_string, ":", simplify = TRUE)
    pair <- str_trim(pair)

    if (length(pair) != 2 || any(!nzchar(pair))) {
      stop(
        "Invalid comparison pair: '", pair_string,
        "'. Use score_type_1:score_type_2."
      )
    }

    tibble(
      score_type_1 = pair[1],
      score_type_2 = pair[2],
      comparison_name = paste0(pair[1], "_minus_", pair[2])
    )
  })

  unknown_types <- setdiff(
    unique(c(parsed$score_type_1, parsed$score_type_2)),
    valid_score_types
  )

  if (length(unknown_types) > 0) {
    stop(
      "Comparison score types were not included in score_types: ",
      paste(unknown_types, collapse = ", ")
    )
  }

  distinct(parsed)
}

comparisons <- parse_comparisons(comparison_arg, score_types)

########################################################
############ Discover phenotypes and scores ############
########################################################

pheno_suffix <- paste0("_", ancestry, ".pheno")

pheno_files <- list.files(
  pheno_dir,
  pattern = paste0("_", ancestry, "\\.pheno$"),
  full.names = TRUE
)

if (length(pheno_files) == 0) {
  stop(
    "No ancestry-specific phenotype files ending in '",
    pheno_suffix,
    "' were found in: ",
    pheno_dir
  )
}

pheno_table <- tibble(
  phenotype = basename(pheno_files) |>
    str_remove(fixed(pheno_suffix)),
  pheno_path = pheno_files
)

parse_score_phenotype <- function(filename, cohort, ancestry) {
  filename |>
    basename() |>
    str_remove(
      paste0("^PG[IS]_", fixed(cohort), "_", fixed(ancestry), "_")
    ) |>
    str_remove("-.*$")
}

score_file_tables <- setNames(
  map(score_types, function(score_type) {
    score_dir <- file.path(pgi_dir, score_type)

    if (!dir.exists(score_dir)) {
      stop("Score-type directory not found: ", score_dir)
    }

    score_prefix_pattern <- paste0(
      "^PG[IS]_",
      cohort,
      "_",
      ancestry,
      "_"
    )

    files <- list.files(
      score_dir,
      pattern = score_prefix_pattern,
      full.names = TRUE,
      recursive = FALSE
    )

    files <- files[!str_detect(files, fixed(paste0(.Platform$file.sep, "parental", .Platform$file.sep)))]

    tibble(
      score_type = score_type,
      phenotype = map_chr(
        files,
        parse_score_phenotype,
        cohort = cohort,
        ancestry = ancestry
      ),
      score_path = files
    ) |>
      group_by(score_type, phenotype) |>
      mutate(number_of_matching_files = n()) |>
      ungroup()
  }),
  score_types
)

duplicate_score_files <- bind_rows(score_file_tables) |>
  filter(number_of_matching_files > 1)

if (nrow(duplicate_score_files) > 0) {
  stop(
    "Multiple score files were found for the same score type and phenotype:\n",
    paste(
      unique(
        paste(
          duplicate_score_files$score_type,
          duplicate_score_files$phenotype,
          sep = " / "
        )
      ),
      collapse = "\n"
    )
  )
}

score_availability <- bind_rows(score_file_tables) |>
  select(score_type, phenotype, score_path) |>
  distinct()

# To preserve a common analysis sample across all requested score types,
# retain phenotypes for which every requested score type exists.
available_counts <- score_availability |>
  count(phenotype, name = "number_of_score_types")

phenotypes_to_predict <- pheno_table |>
  inner_join(available_counts, by = "phenotype") |>
  filter(number_of_score_types == length(score_types)) |>
  arrange(phenotype)

if (nrow(phenotypes_to_predict) == 0) {
  stop(
    "No phenotypes have both a .pheno file and all requested score types."
  )
}

########################################################
#################### Auxiliary data ####################
########################################################

if (cohort == "HRS") {
  if (is.na(crosswalk)) {
    stop("An HRS phenotype-genotype crosswalk must be supplied.")
  }

  score_pheno_crosswalk_data <- read.csv(crosswalk) |>
    transmute(
      IID = as.numeric(SUBJID),
      HHID = as.numeric(HHID),
      PN = as.numeric(PN)
    )
}

if (!is.na(pc_path_arg)) {
  PCs_path <- pc_path_arg
} else if (cohort == "UKB3") {
  PCs_path <- file.path(pheno_dir, "PC_BATCHdum.txt")
} else {
  PCs_path <- file.path(
    pheno_dir,
    "..",
    "..",
    "pc",
    paste0(cohort, "_PCs.eigenvec")
  )
}

if (!file.exists(PCs_path)) {
  stop("PC file not found: ", PCs_path)
}

if (cohort == "UKB3") {
  PCs_data <- fread(PCs_path)
  pc_names <- paste0("PC", 1:20)
  batch_names <- paste0("batch", 1:106)
  covariate_names <- c(pc_names, batch_names)
} else {
  PCs_data <- fread(PCs_path)

  if (ncol(PCs_data) < 22) {
    stop("Expected at least 22 columns in eigenvector file: ", PCs_path)
  }

  old_pc_names <- paste0("V", 3:22)
  new_pc_names <- paste0("pc", 1:20)

  PCs_data <- PCs_data |>
    rename(IID = V2) |>
    rename_with(
      .fn = ~ new_pc_names,
      .cols = all_of(old_pc_names)
    )

  covariate_names <- paste0("pc", 1:10)
}

missing_covariates <- setdiff(covariate_names, names(PCs_data))
if (length(missing_covariates) > 0) {
  stop(
    "Missing expected PC/batch variables: ",
    paste(missing_covariates, collapse = ", ")
  )
}

########################################################
#################### Helper functions ##################
########################################################

read_score_file <- function(score_path, score_type) {
  score_data <- fread(score_path)

  if (!"IID" %in% names(score_data)) {
    stop("Score file lacks an IID column: ", score_path)
  }

  score_column <- if ("SCORE" %in% names(score_data)) {
    "SCORE"
  } else if ("score" %in% names(score_data)) {
    "score"
  } else if (ncol(score_data) >= 5) {
    names(score_data)[5]
  } else {
    stop(
      "Could not identify the score column in: ", score_path,
      ". Expected SCORE, score, or at least five columns."
    )
  }

  score_name <- paste0("score_", score_type)

  score_data |>
    transmute(
      IID = IID,
      !!score_name := as.numeric(.data[[score_column]])
    )
}

prepare_analysis_data <- function(
  cohort,
  pheno_data,
  score_data,
  PCs_data,
  crosswalk_data = NULL
) {
  if (cohort == "HRS") {
    if ("hhidpn" %in% names(pheno_data)) {
      crosswalk_data <- crosswalk_data |>
        mutate(hhidpn = 1000 * HHID + PN)

      analysis_data <- score_data |>
        inner_join(crosswalk_data, by = "IID") |>
        inner_join(pheno_data, by = "hhidpn") |>
        inner_join(PCs_data, by = "IID")
    } else {
      if (
        all(c("hhid", "pn") %in% names(pheno_data)) &&
        !all(c("HHID", "PN") %in% names(pheno_data))
      ) {
        pheno_data <- pheno_data |>
          rename(HHID = hhid, PN = pn)
      }

      if (!all(c("HHID", "PN") %in% names(pheno_data))) {
        stop(
          "HRS phenotype file must contain HHID and PN, ",
          "or legacy hhidpn."
        )
      }

      analysis_data <- score_data |>
        inner_join(crosswalk_data, by = "IID") |>
        inner_join(pheno_data, by = c("HHID", "PN")) |>
        inner_join(PCs_data, by = "IID")
    }
  } else if (cohort == "WLS") {
    if (!"id" %in% names(pheno_data)) {
      stop("WLS phenotype file lacks an id column.")
    }

    analysis_data <- pheno_data |>
      mutate(IID = paste0(id, "_", id)) |>
      inner_join(score_data, by = "IID") |>
      inner_join(PCs_data, by = "IID")
  } else {
    if (!"IID" %in% names(pheno_data)) {
      stop(cohort, " phenotype file lacks an IID column.")
    }

    analysis_data <- pheno_data |>
      inner_join(score_data, by = "IID") |>
      inner_join(PCs_data, by = "IID")
  }

  analysis_data
}

incremental_r2 <- function(data, covariate_names, score_name) {
  baseline_formula <- reformulate(
    termlabels = covariate_names,
    response = "phenotype"
  )

  score_formula <- reformulate(
    termlabels = c(score_name, covariate_names),
    response = "phenotype"
  )

  baseline_model <- lm(baseline_formula, data = data)
  score_model <- lm(score_formula, data = data)

  100 * (
    summary(score_model)$r.squared -
      summary(baseline_model)$r.squared
  )
}

########################################################
################### Run predictions ####################
########################################################

result_rows <- vector("list", nrow(phenotypes_to_predict))
bootstrap_rows <- vector("list", nrow(phenotypes_to_predict))

for (phenotype_index in seq_len(nrow(phenotypes_to_predict))) {
  pheno <- phenotypes_to_predict$phenotype[phenotype_index]
  pheno_path <- phenotypes_to_predict$pheno_path[phenotype_index]

  message("Processing phenotype: ", pheno)

  pheno_data <- fread(pheno_path)

  phenotype_score_files <- score_availability |>
    filter(phenotype == pheno) |>
    arrange(match(score_type, score_types))

  score_data_list <- setNames(
    map(score_types, function(score_type) {
      score_path <- phenotype_score_files |>
        filter(.data$score_type == score_type) |>
        pull(score_path)

      if (length(score_path) != 1) {
        stop(
          "Expected exactly one score file for ", pheno,
          " and ", score_type, "."
        )
      }

      read_score_file(score_path, score_type)
    }),
    score_types
  )

  # Reduce replaces repeated assign()/eval()/parse() calls:
  # it joins every data frame stored in the named score_data_list.
  score_data <- reduce(score_data_list, inner_join, by = "IID")

  analysis_data <- prepare_analysis_data(
    cohort = cohort,
    pheno_data = pheno_data,
    score_data = score_data,
    PCs_data = PCs_data,
    crosswalk_data = if (cohort == "HRS") {
      score_pheno_crosswalk_data
    } else {
      NULL
    }
  )

  score_names <- paste0("score_", score_types)
  required_model_variables <- c(
    "phenotype",
    covariate_names,
    score_names
  )

  missing_model_variables <- setdiff(
    required_model_variables,
    names(analysis_data)
  )

  if (length(missing_model_variables) > 0) {
    stop(
      "Missing model variables for ", pheno, ": ",
      paste(missing_model_variables, collapse = ", ")
    )
  }

  # Fit every baseline and score model in exactly the same complete-case sample.
  analysis_data <- analysis_data |>
    filter(if_all(all_of(required_model_variables), ~ !is.na(.x))) |>
    filter(if_all(all_of(required_model_variables), is.finite))

  if (nrow(analysis_data) == 0) {
    warning("No complete observations for phenotype: ", pheno)
    next
  }

  baseline_formula <- reformulate(
    termlabels = covariate_names,
    response = "phenotype"
  )

  baseline_model <- lm(baseline_formula, data = analysis_data)
  baseline_r2 <- 100 * summary(baseline_model)$r.squared

  observed_incremental_r2 <- setNames(
    map_dbl(score_types, function(score_type) {
      incremental_r2(
        data = analysis_data,
        covariate_names = covariate_names,
        score_name = paste0("score_", score_type)
      )
    }),
    score_types
  )

  bootstrap_incremental_r2 <- matrix(
    NA_real_,
    nrow = iterations,
    ncol = length(score_types),
    dimnames = list(NULL, score_types)
  )

  for (bootstrap_iteration in seq_len(iterations)) {
    set.seed(bootstrap_iteration)

    if (bootstrap_iteration %% 100 == 0) {
      message(
        "Iteration: ", bootstrap_iteration,
        ". Phenotype: ", pheno, "."
      )
    }

    sampled_rows <- sample.int(
      n = nrow(analysis_data),
      size = nrow(analysis_data),
      replace = TRUE
    )
    resampled_data <- analysis_data[sampled_rows, , drop = FALSE]

    bootstrap_incremental_r2[bootstrap_iteration, ] <- map_dbl(
      score_types,
      function(score_type) {
        incremental_r2(
          data = resampled_data,
          covariate_names = covariate_names,
          score_name = paste0("score_", score_type)
        )
      }
    )
  }

  score_summary <- tibble(
    score_type = score_types,
    observed = as.numeric(observed_incremental_r2),
    lower = apply(
      bootstrap_incremental_r2,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper = apply(
      bootstrap_incremental_r2,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    ),
    bootstrap_mean = colMeans(
      bootstrap_incremental_r2,
      na.rm = TRUE
    )
  )

  result_row <- tibble(
    cohort = cohort,
    ancestry = ancestry,
    phenotype = pheno,
    r2_no_score = baseline_r2,
    N = nrow(analysis_data)
  )

  for (score_index in seq_len(nrow(score_summary))) {
    score_type <- score_summary$score_type[score_index]

    result_row[[paste0("r2_inc_", score_type)]] <-
      score_summary$observed[score_index]
    result_row[[paste0("r2_inc_", score_type, "_lower")]] <-
      score_summary$lower[score_index]
    result_row[[paste0("r2_inc_", score_type, "_upper")]] <-
      score_summary$upper[score_index]

    message(
      "Incremental R2 ", score_type, " (95% CI): ",
      round(score_summary$bootstrap_mean[score_index], 4),
      " [",
      round(score_summary$lower[score_index], 4),
      ", ",
      round(score_summary$upper[score_index], 4),
      "]"
    )
  }

  bootstrap_data <- as_tibble(bootstrap_incremental_r2) |>
    rename_with(~ paste0("inc_r2_", .x)) |>
    mutate(
      cohort = cohort,
      ancestry = ancestry,
      iteration = seq_len(iterations),
      phenotype = pheno,
      .before = 1
    )

  if (nrow(comparisons) > 0) {
    for (comparison_index in seq_len(nrow(comparisons))) {
      score_type_1 <- comparisons$score_type_1[comparison_index]
      score_type_2 <- comparisons$score_type_2[comparison_index]
      comparison_name <- comparisons$comparison_name[comparison_index]

      observed_difference <-
        observed_incremental_r2[[score_type_1]] -
        observed_incremental_r2[[score_type_2]]

      bootstrap_difference <-
        bootstrap_incremental_r2[, score_type_1] -
        bootstrap_incremental_r2[, score_type_2]

      comparison_lower <- quantile(
        bootstrap_difference,
        0.025,
        na.rm = TRUE
      )[[1]]
      comparison_upper <- quantile(
        bootstrap_difference,
        0.975,
        na.rm = TRUE
      )[[1]]

      result_row[[paste0("r2_diff_", comparison_name)]] <-
        observed_difference
      result_row[[paste0("r2_diff_", comparison_name, "_lower")]] <-
        comparison_lower
      result_row[[paste0("r2_diff_", comparison_name, "_upper")]] <-
        comparison_upper

      bootstrap_data[[paste0("diff_r2_", comparison_name)]] <-
        bootstrap_difference

      message(
        "Difference ", comparison_name, " (95% CI): ",
        round(mean(bootstrap_difference, na.rm = TRUE), 4),
        " [",
        round(comparison_lower, 4),
        ", ",
        round(comparison_upper, 4),
        "]"
      )
    }
  }

  message("N: ", nrow(analysis_data))

  result_rows[[phenotype_index]] <- result_row
  bootstrap_rows[[phenotype_index]] <- bootstrap_data
}

########################################################
#################### Save results ######################
########################################################

results <- bind_rows(result_rows)
bootstrap_results <- bind_rows(bootstrap_rows)

if (nrow(results) == 0) {
  stop("No phenotype results were produced.")
}

fwrite(
  bootstrap_results,
  file.path(
    bootstrap_output_dir,
    paste0(cohort, "_", ancestry, "_bootstraps.txt")
  ),
  sep = "\t"
)

fwrite(
  results,
  file.path(
    output_dir,
    paste0(cohort, "_", ancestry, "_r2.txt")
  ),
  sep = "\t"
)
