#!/usr/bin/env Rscript

########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "sjmisc", "stringr", "sjmisc", "forcats", "tidyverse", "haven", "magrittr")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)


########################################################
####################### Load data ######################
########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L || length(args) > 4L) {
  stop(
    "Usage: Rscript 10.1.1_construct_WLS_phenotypes.R ",
    "<WLS_renamed.csv> <ancestry> <output_directory> ",
    "[ancestry_FID_IID_file]"
  )
}

input_data_path <- normalizePath(args[[1]], mustWork = TRUE)
ancestry <- toupper(stringr::str_trim(args[[2]]))
output_dir <- normalizePath(args[[3]], mustWork = FALSE)

if (!nzchar(ancestry)) {
  stop("ancestry must be a non-empty label, e.g. EUR, AFR, EAS, SAS, or AMR.")
}

ancestry_fid_iid_path <- if (length(args) == 4L && nzchar(args[[4]])) {
  normalizePath(args[[4]], mustWork = TRUE)
} else {
  file.path(
    dirname(input_data_path),
    paste0("WLS_", ancestry, "_FID_IID.txt")
  )
}

if (!file.exists(input_data_path)) {
  stop("Input phenotype file does not exist: ", input_data_path)
}

if (!file.exists(ancestry_fid_iid_path)) {
  stop("Ancestry ID file does not exist: ", ancestry_fid_iid_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- data.table::fread(input_data_path)

ancestry_ids <- readr::read_delim(
  ancestry_fid_iid_path,
  delim = "\t",
  col_names = c("FID_IID", "IID"),
  col_types = readr::cols(
    FID_IID = readr::col_character(),
    IID = readr::col_character()
  ),
  show_col_types = FALSE
) %>%
  mutate(id = stringr::str_extract(FID_IID, "^[0-9]+")) %>%
  filter(!is.na(id), !stringr::str_detect(id, "[A-Za-z]")) %>%
  pull(id) %>%
  unique()

message("Loaded ", dplyr::n_distinct(data$id), " unique IDs.")
message("Loaded ", length(ancestry_ids), " ", ancestry, " ancestry IDs.")


########################################################
################ Residualising function ################
########################################################

summarise_pheno <- function(data, label = deparse(substitute(data)), preview_n = 6L) {
  if (!all(c("id", "pheno") %in% names(data))) {
    stop(label, " must contain columns named 'id' and 'pheno'.")
  }

  valid_pheno <- data$pheno[!is.na(data$pheno)]

  if (length(valid_pheno) == 0L) {
    summary_stats <- tibble::tibble(
      Min = NA_real_,
      Max = NA_real_,
      Mean = NA_real_,
      Median = NA_real_,
      SD = NA_real_,
      N = dplyr::n_distinct(data$id, na.rm = TRUE)
    )
  } else {
    summary_stats <- data %>%
      summarise(
        Min = min(pheno, na.rm = TRUE),
        Max = max(pheno, na.rm = TRUE),
        Mean = mean(pheno, na.rm = TRUE),
        Median = median(pheno, na.rm = TRUE),
        SD = sd(pheno, na.rm = TRUE),
        N = n_distinct(id, na.rm = TRUE)
      )
  }

  message("\nSummary statistics for ", label, ":")
  print(summary_stats)
  message("First ", preview_n, " rows of ", label, ":")
  print(utils::head(data, preview_n))

  invisible(summary_stats)
}


residualise <- function(data, age_residualise = TRUE, nosex = FALSE) {
  model_formula <- if (age_residualise) {
    pheno ~ age + age2 + male + male_age + male_age2
  } else if (!nosex) {
    pheno ~ dob + dob2 + male + male_dob + male_dob2
  } else {
    pheno ~ dob + dob2
  }

  fit <- stats::lm(model_formula, data = data, na.action = stats::na.exclude)

  data$resid <- stats::residuals(fit)
  data$missing <- as.integer(is.na(data$pheno))

  residual_sd <- stats::sd(data$resid, na.rm = TRUE)
  residual_mean <- mean(data$resid, na.rm = TRUE)

  if (!is.finite(residual_sd) || residual_sd == 0) {
    stop("Residuals cannot be standardised because their SD is zero or undefined.")
  }

  data$std_resid <- (data$resid - residual_mean) / residual_sd
  data
}


residualise.average.save <- function(
  data,
  output_dir,
  ancestry,
  average = TRUE,
  age_residualise = TRUE,
  name,
  nosex = FALSE,
  neb = FALSE
) {
  # Apply the ancestry restriction centrally so every ancestry-labelled
  # phenotype file contains only respondents from the requested ancestry.
  data <- data %>%
    filter(as.character(id) %in% ancestry_ids)

  if (nrow(data) == 0L) {
    stop("No ", ancestry, " ancestry observations remain for phenotype ", name, ".")
  }

  if (average) {
    if (!"wave" %in% names(data)) {
      stop("A 'wave' column is required when average = TRUE for ", name, ".")
    }

    waves <- sort(unique(data$wave[!is.na(data$wave)]))

    df <- lapply(waves, function(current_wave) {
      data_wave <- filter(data, wave == current_wave)
      residualise(
        data = data_wave,
        age_residualise = age_residualise,
        nosex = nosex
      )
    }) %>%
      bind_rows() %>%
      group_by(id, respondent_type) %>%
      summarise(
        phenotype = mean(std_resid, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(is.finite(phenotype))
  } else {
    df <- data %>%
      residualise(
        age_residualise = age_residualise,
        nosex = nosex
      ) %>%
      select(id, respondent_type, phenotype = std_resid)

    if (neb) {
      nebmen <- data %>%
        filter(male == 1) %>%
        residualise(age_residualise = age_residualise, nosex = TRUE) %>%
        select(id, respondent_type, phenotype = std_resid)

      nebwomen <- data %>%
        filter(male == 0) %>%
        residualise(age_residualise = age_residualise, nosex = TRUE) %>%
        select(id, respondent_type, phenotype = std_resid)

      data.table::fwrite(
        nebmen,
        file.path(output_dir, paste0("NEBmen_", ancestry, ".pheno"))
      )
      data.table::fwrite(
        nebwomen,
        file.path(output_dir, paste0("NEBwomen_", ancestry, ".pheno"))
      )
    }
  }

  data.table::fwrite(
    df,
    file.path(output_dir, paste0(name, "_", ancestry, ".pheno"))
  )

  invisible(df)
}


########################################################
############# Construct phenotype: activity ############
########################################################

activity <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("activity")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("activity_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = activity) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(activity)

# residualise, average, save
residualise.average.save(data=activity, average=TRUE, age_residualise=TRUE, name="ACTIVITY", output_dir = output_dir, ancestry = ancestry)


########################################################
############### Construct phenotype: ADHD ##############
########################################################
#dropped for now as ADHD is unavailable for respondents; rather for their children 
# ADHD <- data %>%
#   select(id_old,
#          id,
#          respondent_type,
#          yob,
#          contains("age"),
#          male,
#          african_american,
#          contains("ADHD")) %>%
#   select(-contains("1957"),
#          -contains("1975"),
#          -contains("1993"),
#          -contains("2004")) %>%
#   gather(key = "wave", value = "value",
#                                paste0("ADHD_", 2011),
#                                paste0("age_", 2011)) %>%
#   separate("wave", c("var", "wave")) %>%
#   spread("var", "value") %>%
#   mutate(age2 = age^2,
#          male_age = male * age,
#          male_age2 = male * age2,
#          pheno = ADHD) %>% # reverse coded in construction stage
#   select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
#   drop_na()

# residualise, average, save
# residualise.average.save(data=ADHD, average=FALSE, age_residualise=TRUE, name="ADHD", output_dir = output_dir, ancestry = ancestry)


########################################################
######### Construct phenotype: age first birth #########
########################################################

AFB <- data %>%
  select(
    id_old,
    id,
    respondent_type,
    yob,
    contains("age"),
    male,
    african_american,
    contains("AFB")
  ) %>%
  select(-contains("1957"), -contains("1993"), -contains("2004")) %>%
  mutate(
    AFB = case_when(
      !is.na(AFB_2011) & AFB_2011 > 10 ~ AFB_2011,  # Only accept AFB_2011 values greater than 10
      AFB_1975 > 10 ~ AFB_1975,                    # Only accept AFB_1975 values greater than 10
      TRUE ~ NA_real_                             # Set to NA if neither condition is met
    ),
    dob = yob,
    dob2 = yob^2,
    male_dob = male * dob,
    male_dob2 = male * dob2,
    pheno = AFB  # Use refined AFB values for pheno
  ) %>%
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
AFB <- AFB %>%
  filter(id %in% ancestry_ids)

summarise_pheno(AFB)

# Add checks to understand the values that go into pheno calculation
check_values <- data %>%
  summarise(
    Max_AFB_2011 = max(AFB_2011, na.rm = TRUE),
    Min_AFB_2011 = min(AFB_2011, na.rm = TRUE),
    Max_AFB_1975 = max(AFB_1975, na.rm = TRUE),
    Min_AFB_1975 = min(AFB_1975, na.rm = TRUE)
  )

# Print the check values
print("Check values for AFB_2011 and AFB_1975:")
print(check_values)

# residualise, average, save
residualise.average.save(data=AFB, average=FALSE, age_residualise=FALSE, name="AFB", output_dir = output_dir, ancestry = ancestry)

########################################################
#### Construct phenotype: age at smoking initiation ####
########################################################

ASI <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         african_american,
         contains("startSmoke")) %>%
  gather(key="wave", value="value",
         paste0("startSmoke_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(ASI_average = mean(startSmoke, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob   = yob,
         dob2  = yob^2,
         pheno = ASI_average) %>% # reverse code
  select(id, respondent_type, wave, pheno, dob, dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
ASI <- ASI %>%
  filter(id %in% ancestry_ids)

summarise_pheno(ASI)

# residualise, average, save
residualise.average.save(data=ASI, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="ASI", output_dir = output_dir, ancestry = ancestry)

########################################################
######### Construct phenotype: age first menses ########
########################################################

AFM <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         african_american,
         contains("AFM")) %>%
  select(-contains("1957"), -contains("1975"), -contains("2011")) %>%
  gather(key="wave", value="value",
         paste0("AFM_", c(1993, 2004)), paste0("age_", c(1993, 2004))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(AFM_average = mean(AFM, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob   = yob,
         dob2  = yob^2,
         pheno = AFM_average) %>% # reverse code
  select(id, respondent_type, wave, pheno, dob, dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
AFM <- AFM %>%
  filter(id %in% ancestry_ids)

summarise_pheno(AFM)

# residualise, average, save
residualise.average.save(data=AFM, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="MENARCHE", output_dir = output_dir, ancestry = ancestry)

########################################################
########## Construct phenotype: agreeableness ##########
########################################################

#agree <- data %>%
#  select(id_old,
#         id,
#         respondent_type,
#         yob,
#         contains("age"),
#         male,
#         african_american,
#         contains("agree")) %>%
#  select(-contains("phone"), -contains("nanswered")) %>%
#  select(-contains("1957"), -contains("1975")) %>%
#  gather(key = "wave", value = "value",
#                               paste0("agree_", c(1993, 2004, 2011)),
#                               paste0("age_",   c(1993, 2004, 2011))) %>%
#  separate("wave", c("var", "wave")) %>%
#  spread("var", "value") %>%
#  mutate(age2 = age^2,
#         male_age = male * age,
#         male_age2 = male * age2,
#         pheno = agree) %>%
#  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
#  drop_na()

# residualise, average, save
#residualise.average.save(data=agree, average=TRUE, age_residualise=TRUE, name="AGREE", output_dir = output_dir, ancestry = ancestry)

########################################################
############## Construct phenotype: Alzheimer's #############
########################################################

alzheimers <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("alzheimers")) %>%
  gather(key="wave", value="value",
         paste0("alzheimers_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(alzheimers = max(alzheimers, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = alzheimers) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
alzheimers <- alzheimers %>%
  filter(id %in% ancestry_ids)

summarise_pheno(alzheimers)

# residualise, average, save
residualise.average.save(data=alzheimers, average=FALSE, age_residualise=FALSE, name="ALZ", output_dir = output_dir, ancestry = ancestry)

########################################################
############## Construct phenotype: asthma #############
########################################################

asthma <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("asthma"),
         -contains("hayfever")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("asthma_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  filter(!all(is.na(asthma))) %>%  # Exclude groups where all asthma values are NA
  mutate(asthma = max(asthma, na.rm=TRUE),
         row_n = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = asthma) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
asthma <- asthma %>%
  filter(id %in% ancestry_ids)

summarise_pheno(asthma)

# residualise, average, save
residualise.average.save(data=asthma, average=FALSE, age_residualise=FALSE, name="ASTHMA", output_dir = output_dir, ancestry = ancestry)

########################################################
########## Construct phenotype: asthmahayfever #########
########################################################

asthmahayfever <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("asthma_2011"),
         contains("hayfever_2011")) %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(asthma_2011 + hayfever_2011 >  0 ~ 1,
                           asthma_2011 + hayfever_2011 == 0 ~ 0)) %>%
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
asthmahayfever <- asthmahayfever %>%
  filter(id %in% ancestry_ids)

summarise_pheno(asthmahayfever)

# residualise, average, save
residualise.average.save(data=asthmahayfever, average=FALSE, age_residualise=FALSE, name="ASTECZRHI", output_dir = output_dir, ancestry = ancestry)


########################################################
############## Construct phenotype: audit ##############
########################################################

audit <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("audit")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("audit_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = audit) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(audit)

# residualise, average, save
residualise.average.save(data=audit, average=T, age_residualise=T, name="AUDIT", output_dir = output_dir, ancestry = ancestry)


########################################################
############### Construct phenotype: bmi ###############
########################################################

# Step 1: Gather and clean BMI data
bmi <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("bmi")) %>%
  select(-contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("bmi_", c(1957, 1993, 2004, 2011)), paste0("age_", c(1957, 1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = bmi) %>%
  filter(pheno >= 10) %>% # biologically implausible to have bmi <10
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Step 2: Calculate Z-scores and filter based on CDC guidelines
# Calculate mean and standard deviation for z-score computation
mean_bmi <- mean(bmi$pheno, na.rm = TRUE)
sd_bmi <- sd(bmi$pheno, na.rm = TRUE)

# Adding a temporary z-score column to filter out biologically implausible values
bmi <- bmi %>%
  mutate(z_score = (pheno - mean_bmi) / sd_bmi) %>%
  filter(z_score >= -6, z_score <= 6) %>%
  select(-z_score)  # Remove the z-score column after filtering

summarise_pheno(bmi)

# residualise, average, save
residualise.average.save(data=bmi, average=TRUE, age_residualise=TRUE, name="BMI", output_dir = output_dir, ancestry = ancestry)

########################################################
######### Construct phenotype: Breast cancer ###########
########################################################

BRCA <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("breastCancer")) %>%
  gather(key="wave", value="value",
         paste0("breastCancer_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(breastCancer = max(breastCancer, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  filter(male == 0) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = breastCancer) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
BRCA <- BRCA %>%
  filter(id %in% ancestry_ids)

summarise_pheno(BRCA)

# residualise, average, save
residualise.average.save(data=BRCA, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BRCA", output_dir = output_dir, ancestry = ancestry)

########################################################
########### Construct phenotype: cat allergy ###########
########################################################

cat <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("cat")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993"), -contains("2004")) %>%
  gather(key="wave", value="value",
         paste0("cat_", 2011), paste0("age_", 2011)) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = cat) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
cat <- cat %>%
  filter(id %in% ancestry_ids)

summarise_pheno(cat)

# residualise, average, save
residualise.average.save(data=cat, average=FALSE, age_residualise=TRUE, name="ALLERGYCAT", output_dir = output_dir, ancestry = ancestry)


########################################################
########## Construct phenotype: conscientious ##########
########################################################

consc <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("consc"), -contains("phone")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("consc_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = consc) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(consc)

# residualise, average, save
residualise.average.save(data=consc, average=TRUE, age_residualise=TRUE, name="CONSC", output_dir = output_dir, ancestry = ancestry)


########################################################
############### Construct phenotype: COPD ##############
########################################################

copd <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("copd"), -contains("phone")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("copd_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(copd = ifelse(any(!is.na(copd)), max(copd, na.rm=TRUE), NA),
         row_n = row_number()) %>%
  ungroup() %>%
  filter(row_n == 1 & copd >= 0) %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = copd) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
copd <- copd %>%
  filter(id %in% ancestry_ids)

summarise_pheno(copd)

# residualise, average, save
residualise.average.save(data=copd, average=FALSE, age_residualise=FALSE, name="COPD", output_dir = output_dir, ancestry = ancestry)

########################################################
########### Construct phenotype: cigs per day ##########
########################################################

CPD <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("CPD"), -contains("phone")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("CPD_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = CPD) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
CPD <- CPD %>%
  filter(id %in% ancestry_ids)

summarise_pheno(CPD)

# residualise, average, save
residualise.average.save(data=CPD, average=TRUE, age_residualise=TRUE, name="CPD", output_dir = output_dir, ancestry = ancestry)


########################################################
############ Construct phenotype: depression ###########
########################################################

depr <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("depr"), -contains("nanswered")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("depr_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = depr) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(depr)

# residualise, average, save
residualise.average.save(data=depr, average=TRUE, age_residualise=TRUE, name="DEP", output_dir = output_dir, ancestry = ancestry)

########################################################
############### Construct phenotype: DPW ###############
########################################################

dpw <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("dpw"), -contains("nanswered")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("dpw_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = dpw) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(dpw)

# residualise, average, save
residualise.average.save(data=dpw, average=TRUE, age_residualise=TRUE, name="DPW", output_dir = output_dir, ancestry = ancestry)


########################################################
########### Construct phenotype: dust allergy ##########
########################################################

dust <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("dust")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993"), -contains("2004")) %>%
  gather(key="wave", value="value",
         paste0("dust_", 2011), paste0("age_", 2011)) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = dust) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
dust <- dust %>%
  filter(id %in% ancestry_ids)

summarise_pheno(dust)

# residualise, average, save
residualise.average.save(data=dust, average=FALSE, age_residualise=TRUE, name="ALLERGYDUST", output_dir = output_dir, ancestry = ancestry)


########################################################
################ Construct phenotype: EA ###############
########################################################

EA <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("EA_")) %>%
  #mutate(EA_diff = EA_2011 - EA_2004) %>% # EA values barely change, as expected
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(!is.na(EA_2011) ~ EA_2011, TRUE ~ EA_2004)) %>%
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
EA <- EA %>%
  filter(id %in% ancestry_ids)

summarise_pheno(EA)

# residualise, average, save
residualise.average.save(data=EA, average=FALSE, age_residualise=FALSE, name="EA", output_dir = output_dir, ancestry = ancestry)


########################################################
########### Construct phenotype: ever smoker ###########
########################################################

eversmoke <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("eversmoke")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("eversmoke_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  filter(!all(is.na(eversmoke))) %>%  # Exclude groups where all eversmoke values are NA
  mutate(eversmoke_max = max(eversmoke, na.rm=TRUE),
         row_n = row_number()) %>%
  ungroup() %>%
  filter(row_n == 1 & eversmoke_max >= 0) %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = eversmoke_max) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
eversmoke <- eversmoke %>%
  filter(id %in% ancestry_ids)

summarise_pheno(eversmoke)

# residualise, average, save
residualise.average.save(data=eversmoke, average=FALSE, age_residualise=FALSE, name="EVERSMOKE", output_dir = output_dir, ancestry = ancestry)


########################################################
########### Construct phenotype: extraversion ##########
########################################################

extra <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("extra"), -contains("nanswered"), -contains("phone")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("extra_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = extra) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(extra)

# residualise, average, save
residualise.average.save(data=extra, average=TRUE, age_residualise=TRUE, name="EXTRA", output_dir = output_dir, ancestry = ancestry)


########################################################
####### Construct phenotype: family satisfaction #######
########################################################

famsat <- data %>%
  select(id_old, id, respondent_type, yob, contains("age_2004"), male, african_american, contains("famsat")) %>%
  mutate(age = age_2004,
         age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = famsat_2004) %>% # reverse code
  select(id, respondent_type, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(famsat)

# residualise, average, save
residualise.average.save(data=famsat, average=FALSE, age_residualise=TRUE, name="FAMSAT", output_dir = output_dir, ancestry = ancestry)


########################################################
###### Construct phenotype: financial satisfaction #####
########################################################

finsat <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("finsat")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993")) %>%
  gather(key="wave", value="value",
         paste0("finsat_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = finsat) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(finsat)

# residualise, average, save
residualise.average.save(data=finsat, average=TRUE, age_residualise=TRUE, name="FINSAT", output_dir = output_dir, ancestry = ancestry)


########################################################
####### Construct phenotype: friend satisfaction #######
########################################################

friendsat1 <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("friendsat1")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993"), -contains("2011")) %>%
  mutate(age = age_2004, friendsat1 = friendsat1_2004, wave = 2004) %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = friendsat1) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

friendsat2 <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("friendsat2")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993"), -contains("2011")) %>%
  mutate(age = age_2004, friendsat2 = friendsat2_2004, wave = 2004) %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = friendsat2) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(friendsat1)

# residualise, average, save
residualise.average.save(data=friendsat1, average=TRUE, age_residualise=TRUE, name="FRIENDSAT1", output_dir = output_dir, ancestry = ancestry)
residualise.average.save(data=friendsat2, average=TRUE, age_residualise=TRUE, name="FRIENDSAT2", output_dir = output_dir, ancestry = ancestry)


########################################################
############ Construct phenotype: hay-fever ############
########################################################

hayfever <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("hayfever"), -contains("asthma")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993"), -contains("2004")) %>%
  gather(key="wave", value="value",
         paste0("hayfever_", 2011), paste0("age_", 2011)) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = hayfever) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
hayfever <- hayfever %>%
  filter(id %in% ancestry_ids)

summarise_pheno(hayfever)

# residualise, average, save
residualise.average.save(data=hayfever, average=FALSE, age_residualise=FALSE, name="HAYFEVER", output_dir = output_dir, ancestry = ancestry)

###########################################################################################################
##################### Construct phenotype: heart attack or mycardial infraction (for HARDCAD) #############
###########################################################################################################

HEARTATTACK <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("heartattack")) %>%
  gather(key="wave", value="value",
          paste0("heartattack_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  filter(!all(is.na(heartattack))) %>%  # Exclude groups where all heartattack values are NA
  mutate(heartattack = max(heartattack, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = heartattack) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  filter(is.finite(pheno)) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
HEARTATTACK <- HEARTATTACK %>%
  filter(id %in% ancestry_ids)

summarise_pheno(HEARTATTACK)

# residualise, average, save
residualise.average.save(data=HEARTATTACK, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="HARDCAD", output_dir = output_dir, ancestry = ancestry)


########################################################
############## Construct phenotype: height #############
########################################################

height <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("height")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("height_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(height = mean(height, na.rm=TRUE),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn == 1) %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = height) %>%
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

summarise_pheno(height)

# residualise, average, save
residualise.average.save(data=height, average=FALSE, age_residualise=FALSE, name="HEIGHT", output_dir = output_dir, ancestry = ancestry)

##########################################################################
####### Construct phenotype: high cholesterol (used for BL_CHOL, BL_HDL, BL_LDL, BL_TRYG)###########
##########################################################################

HIGHCHOL <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("highcholesterol")) %>%
  gather(key="wave", value="value",
          paste0("highcholesterol_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  filter(!all(is.na(highcholesterol))) %>%  # Exclude groups where all highcholesterol values are NA
  mutate(highcholesterol = max(highcholesterol, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = highcholesterol) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  filter(is.finite(pheno)) %>%
  drop_na()

summarise_pheno(HIGHCHOL)

# residualise, average, save
residualise.average.save(data=HIGHCHOL, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BL_CHOL", output_dir = output_dir, ancestry = ancestry)
residualise.average.save(data=HIGHCHOL, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BL_HDL", output_dir = output_dir, ancestry = ancestry)
residualise.average.save(data=HIGHCHOL, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BL_LDL", output_dir = output_dir, ancestry = ancestry)
residualise.average.save(data=HIGHCHOL, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BL_TRYG", output_dir = output_dir, ancestry = ancestry)

####################################################################################
####### Construct phenotype: high blood pressure (used for BPdia, BPsys)###########
#####################################################################################

highBP <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("highBP")) %>%
  gather(key="wave", value="value",
          paste0("highBP_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  filter(!all(is.na(highBP))) %>%  # Exclude groups where all highcholesterol values are NA
  mutate(highBP = max(highBP, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = highBP) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  filter(is.finite(pheno)) %>%
  drop_na()

summarise_pheno(highBP)

# residualise, average, save
residualise.average.save(data=highBP, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BPdia", output_dir = output_dir, ancestry = ancestry)
residualise.average.save(data=highBP, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="BPsys", output_dir = output_dir, ancestry = ancestry)


########################################################
######### Construct phenotype: Insomina ################
########################################################

INSOMNIA <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("insomnia")) %>%
  gather(key="wave", value="value",
         paste0("insomnia_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(insomnia = max(insomnia, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = insomnia) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

summarise_pheno(INSOMNIA)

# residualise, average, save
residualise.average.save(data=INSOMNIA, average=FALSE, age_residualise=FALSE, name="INSOMNIA", output_dir = output_dir, ancestry = ancestry)

########################################################
########### Construct phenotype: intelligence ##########
########################################################

intelligence <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("intelligence_")) %>%
  gather(key="wave", value="value", paste0("intelligence_", c(1957, 1975)), paste0("age_", c(1957, 1975))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = intelligence) %>%
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

summarise_pheno(intelligence)

# residualise, average, save
residualise.average.save(data=intelligence, average=FALSE, age_residualise=FALSE, name="CP", output_dir = output_dir, ancestry = ancestry)


########################################################
######### Construct phenotype: left out social #########
########################################################

#leftoutsocial <- data %>%
#  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("leftoutsocial"), -contains("phone")) %>%
#  select(-contains("1957"), -contains("1975")) %>%
#  gather(key="wave", value="value",
#         paste0("leftoutsocial_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
#  separate("wave", c("var", "wave")) %>%
#  spread("var", "value") %>%
#  mutate(age2 = age^2,
#         male_age = male * age,
#         male_age2 = male * age2,
#         pheno = leftoutsocial) %>%
#  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
#  drop_na()

# residualise, average, save
#residualise.average.save(data=leftoutsocial, average=TRUE, age_residualise=TRUE, name="LEFTOUT", output_dir = output_dir, ancestry = ancestry)


########################################################
############## Construct phenotype: lonely #############
########################################################

lonely <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("lonely")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("lonely_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = lonely) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(lonely)

# residualise, average, save
residualise.average.save(data=lonely, average=TRUE, age_residualise=TRUE, name="LONELY", output_dir = output_dir, ancestry = ancestry)


########################################################
############# Construct phenotype: migraine ############
########################################################

#migraine <- data %>%
#  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("migraine")) %>%
#  select(-contains("1957"), -contains("1975")) %>%
#  gather(key="wave", value="value",
#         paste0("migraine_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
#  separate("wave", c("var", "wave")) %>%
#  spread("var", "value") %>%
#  mutate(age2 = age^2,
#         male_age = male * age,
#         male_age2 = male * age2,
#         pheno = migraine) %>%
#  filter(pheno >= 0) %>%
#  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
#  drop_na()

# residualise, average, save
#residualise.average.save(data=migraine, average=TRUE, age_residualise=TRUE, name="MIGRAINE", output_dir = output_dir, ancestry = ancestry)


########################################################
######### Construct phenotype: number ever born ########
########################################################

NEB <- data %>%
  select(id_old, id, respondent_type, yob, male, african_american, NEB, contains("age")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key = "wave", value = "age_value", paste0("age_", c("1993", "2004", "2011"))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "age_value") %>%
  mutate(age = as.numeric(age), # Convert the 'age_value' to numeric and rename as 'age'
         dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = NEB) %>%
  select(id, respondent_type, wave, pheno, age, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Calculate minimum age across waves for each individual
NEB <- NEB %>%
  group_by(id) %>%
  mutate(min_age = min(age, na.rm = TRUE)) %>%
  ungroup() %>%
  filter((male == 0 & min_age >= 45) | (male == 1 & min_age >= 55)) # Apply age restrictions based on minimum age

       # Diagnostic for age-based filtering
       age_filter_check <- with(NEB, table(male, Age_Group = ifelse(male == 1, min_age >= 55, min_age >= 45)))
       print(age_filter_check)
       # Additional diagnostics for minimum age
       min_age_female <- min(NEB$min_age[NEB$male == 0], na.rm = TRUE)
       print(paste("Minimum age for females:", min_age_female))
       min_age_male <- min(NEB$min_age[NEB$male == 1], na.rm = TRUE)
       print(paste("Minimum age for males:", min_age_male))

# Filter the dataset to include only the requested ancestry IDs
NEB <- NEB %>%
  filter(id %in% ancestry_ids)

# Restrict by gender to get NEBmen and NEBwomen stats 
#NEB <- NEB %>%
#  filter(id %in% ancestry_ids, male == 1)

summarise_pheno(NEB)

# Residualise, average, save
residualise.average.save(data=NEB, average=FALSE, age_residualise=FALSE, name="NEBpooled", neb=TRUE, output_dir = output_dir, ancestry = ancestry)


########################################################
########### Construct phenotype: neuroticism ###########
########################################################

neur <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("neur")) %>%
  select(-contains("phone"), -contains("nanswered")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("neur_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = neur) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(neur)

# residualise, average, save
residualise.average.save(data=neur, average=TRUE, age_residualise=TRUE, name="NEURO", output_dir = output_dir, ancestry = ancestry)


########################################################
############# Construct phenotype: openness ############
########################################################

open <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("open"), -contains("phone")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("open_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = open) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(open)

# residualise, average, save
residualise.average.save(data=open, average=TRUE, age_residualise=TRUE, name="OPEN", output_dir = output_dir, ancestry = ancestry)


########################################################
######### Construct phenotype: prostate cancer ###########
########################################################

PRCA <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("prostateCancer")) %>%
  gather(key="wave", value="value",
         paste0("prostateCancer_", c(2004, 2011)), paste0("age_", c(2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(prostateCancer = max(prostateCancer, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  filter(male == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = prostateCancer) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
PRCA <- PRCA %>%
  filter(id %in% ancestry_ids)

summarise_pheno(PRCA)

# residualise, average, save
residualise.average.save(data=PRCA, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="PRCA", output_dir = output_dir, ancestry = ancestry)


########################################################
########## Construct phenotype: pollen allergy #########
########################################################

pollen <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("pollen")) %>%
  select(-contains("1957"), -contains("1975"), -contains("1993"), -contains("2004")) %>%
  gather(key="wave", value="value",
         paste0("pollen_", 2011), paste0("age_", 2011)) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = pollen) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
pollen <- pollen %>%
  filter(id %in% ancestry_ids)

summarise_pheno(pollen)

# residualise, average, save
residualise.average.save(data=pollen, average=FALSE, age_residualise=TRUE, name="ALLERGYPOLLEN", output_dir = output_dir, ancestry = ancestry)


########################################################
########### Construct phenotype: religiosity ###########
########################################################

relig <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("relig")) %>%
  select(-contains("1957"), -contains("phone"), -contains("mail")) %>%
  gather(key="wave", value="value",
         paste0("relig_", c(1975, 1993, 2004, 2011)), paste0("age_", c(1975, 1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = relig) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(relig)

# residualise, average, save
residualise.average.save(data=relig, average=TRUE, age_residualise=TRUE, name="RELIGATT", output_dir = output_dir, ancestry = ancestry)


########################################################
############### Construct phenotype: risk ##############
########################################################

risk <- data %>%
  select(id_old, id, respondent_type, yob, contains("age_2011"), male, african_american, contains("risk")) %>%
  select(-contains("losing")) %>%
  gather(key="wave", value="value",
         paste0("risk", c(5, 9, 11), "_2011")) %>%
  separate("wave", c("wave", "var")) %>%
  mutate(age = age_2011,
         age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = value,
         wave = str_remove(wave, "risk")) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

risk_loss <- data %>%
  select(id_old, id, respondent_type, yob, contains("age_2011"), male, african_american, contains("risklosing")) %>%
  gather(key="wave", value="value",
         paste0("risklosing", c(5, 9, 11), "_2011")) %>%
  separate("wave", c("wave", "var")) %>%
  mutate(age = age_2011,
         age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = value,
         wave = str_remove(wave, "risklosing")) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(risk)

# residualise, average, save
residualise.average.save(data=risk, average=TRUE, age_residualise=TRUE, name="RISK", output_dir = output_dir, ancestry = ancestry)
residualise.average.save(data=risk_loss, average=TRUE, age_residualise=TRUE, name="RISKLOSS", output_dir = output_dir, ancestry = ancestry)


########################################################
######## Construct phenotype: self-rated health ########
########################################################

selfhealth <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("selfhealth")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("selfhealth_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = selfhealth) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(selfhealth)

# residualise, average, save
residualise.average.save(data=selfhealth, average=TRUE, age_residualise=TRUE, name="SELFHEALTH", output_dir = output_dir, ancestry = ancestry)


########################################################
######### Construct phenotype: smoking cessation ###########
########################################################
SMCESS <- data %>%
   select(id_old,
          id,
          respondent_type,
          yob,
          contains("age"),
          male,
          african_american,
          contains("formerSmoker")) %>%
   gather(key="wave", value="value",
          paste0("formerSmoker_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
   separate("wave", c("var", "wave")) %>%
   spread("var", "value") %>%
   group_by(id, respondent_type) %>%
   filter(!all(is.na(formerSmoker))) %>%  # Exclude groups where all formerSmoker values are NA
   mutate(formerSmoker = max(formerSmoker, na.rm=TRUE),
          row_n       = row_number()) %>%
   filter(row_n == 1) %>%
   ungroup() %>%
   mutate(dob = yob,
          dob2 = yob^2,
          male_dob = male * dob,
          male_dob2 = male * dob2,
          pheno = formerSmoker) %>%
   filter(pheno >= 0) %>%
   select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
   drop_na()

# Filter the dataset to include only the requested ancestry IDs
SMCESS <- SMCESS %>%
  filter(id %in% ancestry_ids)

summarise_pheno(SMCESS)

# residualise, average, save
residualise.average.save(data=SMCESS, average=FALSE, age_residualise=FALSE, name="SMCESS", output_dir = output_dir, ancestry = ancestry)


########################################################
###### Construct phenotype: subjective well-being ######
########################################################

SWB <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("SWB")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("SWB_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = SWB) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(SWB)

# residualise, average, save
residualise.average.save(data=SWB, average=TRUE, age_residualise=TRUE, name="SWB", output_dir = output_dir, ancestry = ancestry)


########################################################
####### Construct phenotype: type-2 diabetes ###########
########################################################

T2D <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("type2diabetes")) %>%
  gather(key="wave", value="value",
          paste0("type2diabetes_", 2011), paste0("age_", 2011)) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  filter(!all(is.na(type2diabetes))) %>%  # Exclude groups where all type2diabetes values are NA
  mutate(type2diabetes = max(type2diabetes, na.rm=TRUE),
         row_n       = row_number()) %>%
  filter(row_n == 1) %>%
  ungroup() %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2, 
         pheno = type2diabetes) %>%
  select(id, respondent_type, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  filter(is.finite(pheno)) %>%
  drop_na()

# Filter the dataset to include only the requested ancestry IDs
T2D <- T2D %>%
  filter(id %in% ancestry_ids)
  
summarise_pheno(T2D)

# residualise, average, save
residualise.average.save(data=T2D, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="T2D", output_dir = output_dir, ancestry = ancestry)

########################################################
######## Construct phenotype: work satisfaciton ########
########################################################

worksat <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("worksat")) %>%
  select(-contains("1957")) %>%
  gather(key="wave", value="value",
         paste0("worksat_", c(1975, 1993, 2004, 2011)), paste0("age_", c(1975, 1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = worksat) %>% # reverse code
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

summarise_pheno(worksat)

# residualise, average, save
residualise.average.save(data=worksat, average=TRUE, age_residualise=TRUE, name="WORKSAT", output_dir = output_dir, ancestry = ancestry)
