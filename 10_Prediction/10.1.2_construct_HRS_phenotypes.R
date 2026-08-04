#!/usr/bin/env Rscript

########################################################
######################## Set-up ########################
########################################################


# Command-line arguments
#
# 1. inputDataDir: directory containing RAND-HRS, wave, tracker, and
#    ancestry-ID input files
# 2. outputDataDir: directory where phenotype files will be written
# 3. ancestry: ancestry label, e.g. EUR, AFR, EAS, SAS, or AMR
# 4. libraryDir: optional personal R-library directory
# 5. ancestryIdFile: optional path to the ancestry-specific FID/IID file
#
# When ancestryIdFile is omitted, the script expects:
#   <inputDataDir>/HRS_<ancestry>_FID_IID.txt
#
# Example:
# Rscript 10.1.2_construct_HRS_phenotypes.R \
#   /path/to/input \
#   /path/to/output \
#   AFR \
#   /path/to/R/library \
#   /path/to/HRS_AFR_FID_IID.txt

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 5) {
  stop(
    "Usage: Rscript script.R <inputDataDir> <outputDataDir> <ancestry> ",
    "[libraryDir] [ancestryIdFile]"
  )
}

inputDataDir <- normalizePath(args[1], mustWork = TRUE)
outputDataDir <- normalizePath(args[2], mustWork = FALSE)
ancestry <- toupper(trimws(args[3]))

if (!nzchar(ancestry)) {
  stop("ancestry must be a non-empty label, e.g. EUR, AFR, EAS, SAS, or AMR.")
}

libraryDir <- if (length(args) >= 4 && nzchar(args[4]) && toupper(args[4]) != "NA") {
  normalizePath(args[4], mustWork = FALSE)
} else {
  NULL
}

ancestryIdFile <- if (length(args) >= 5 && nzchar(args[5]) && toupper(args[5]) != "NA") {
  normalizePath(args[5], mustWork = TRUE)
} else {
  file.path(inputDataDir, paste0("HRS_", ancestry, "_FID_IID.txt"))
}

if (!dir.exists(outputDataDir)) {
  dir.create(outputDataDir, recursive = TRUE)
}

if (!is.null(libraryDir)) {
  if (!dir.exists(libraryDir)) {
    dir.create(libraryDir, recursive = TRUE)
  }
  .libPaths(c(libraryDir, .libPaths()))
}

# List of packages to install and load.
packages <- c(
  "dplyr", "tidyverse", "data.table", "foreign",
  "sjmisc", "haven", "purrr", "tidyr"
)

for (pkg in unique(packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      lib = if (is.null(libraryDir)) .libPaths()[1] else libraryDir
    )
  }
  library(pkg, character.only = TRUE)
}

########################################################
#################### Load RAND data ####################
########################################################

RAND_path <- file.path(inputDataDir, "randhrs1992_2016v2.dta")
RAND_data <- read.dta(RAND_path)


# Load ancestry-specific household IDs.
if (!file.exists(ancestryIdFile)) {
  stop("Ancestry ID file not found: ", ancestryIdFile)
}

ancestry_hhids <- readr::read_delim(
  ancestryIdFile,
  delim = "\t",
  col_names = c("HHID", "IID"),
  col_types = readr::cols(
    HHID = readr::col_character(),
    IID = readr::col_double()
  ),
  show_col_types = FALSE
) %>%
  mutate(HHID = as.integer(stringr::str_remove(HHID, "^0+"))) %>%
  pull(HHID)


   # Checking for write permission in outputDataDir
    has_write_perm <- file.access(outputDataDir, 2) == 0

    if (has_write_perm) {
      print("You have write permission in this directory.")
    } else {
      print("You do not have write permission in this directory.")
    }

########################################################
################ Residualising functions ###############
########################################################

# Function to perform residualization and standardization of phenotype data.
residualise <- function(data, age_residualise=TRUE, sexres=TRUE) {
  # Decide on the regression model based on the input flags for age and sex adjustments.
  if (age_residualise && sexres) {
    # Model incorporating both age and sex interactions.
    reg <- summary(lm(pheno ~ age + I(age^2) + male + male:age + male:I(age^2), data))
  } else if (age_residualise) {
    # Model including only age and age squared.
    reg <- summary(lm(pheno ~ age + I(age^2), data))
  } else if (sexres) {
    # Model adjusting for demographic variables like dob and sex.
    reg <- summary(lm(pheno ~ dob + I(dob^2) + male + male:dob + male:I(dob^2), data))
  } else {
    # Simplified model with only dob and dob squared.
    reg <- summary(lm(pheno ~ dob + I(dob^2), data))
  }

  # Identify non-missing phenotype observations for residual calculations.
  sel <- which(!is.na(data$pheno))
  data$resid <- NA
  data$resid[sel] <- reg$residuals
  data$missing <- 1
  data$missing[sel] <- 0

  # Standardize the residuals for consistent scale across all data.
  mean_resid <- mean(data$resid, na.rm = TRUE)
  sd_resid <- sd(data$resid, na.rm = TRUE)
  data$std_resid <- (data$resid - mean_resid) / sd_resid

  return(data)
}

# Function to process data, perform residualization, and optionally average and save the results.
residualise.average.save <- function(data, average=TRUE, age_residualise=TRUE, name=NULL, sexres=TRUE, rand=FALSE, save=TRUE) {
  # Apply the ancestry restriction centrally so every ancestry-labelled
  # phenotype file contains only respondents from the requested ancestry.
  if (rand) {
    if (!"hhidpn" %in% names(data)) {
      stop("RAND phenotype data must contain hhidpn for ancestry filtering.")
    }
    data <- data %>%
      mutate(.ancestry_HHID = as.numeric(substr(as.character(hhidpn), 1, 5))) %>%
      filter(.ancestry_HHID %in% ancestry_hhids) %>%
      select(-.ancestry_HHID)
  } else {
    if (!"HHID" %in% names(data)) {
      stop("HRS phenotype data must contain HHID for ancestry filtering.")
    }
    data <- data %>%
      filter(HHID %in% ancestry_hhids)
  }

  if (nrow(data) == 0L) {
    stop("No ", ancestry, " ancestry observations remain for phenotype ", name, ".")
  }

  # Initialize an empty dataframe to hold processed results.
  df <- data.frame()

  # Check if the data should be processed wave by wave.
  if (average) {
    # Identify unique waves in the dataset for separate processing.
    waves <- unique(data$wave)
    for (wave_value in waves) {
      # Filter data by wave and apply the residualise function.
      data_wave <- filter(data, .data$wave == wave_value)
      data_wave <- residualise(data_wave, age_residualise, sexres)
      # Combine processed data from each wave.
      df <- rbind(df, data_wave)
    }

    # Average standardized residuals if specified, either by household ID or individual participant number.
    if (rand) {
      df <- df %>%
        group_by(hhidpn) %>%
        summarise(phenotype = mean(std_resid, na.rm = TRUE)) %>%
        filter(!is.na(phenotype))
    } else {
      df <- df %>%
        group_by(HHID, PN) %>%
        summarise(phenotype = mean(std_resid, na.rm = TRUE)) %>%
        filter(!is.na(phenotype))
    }
  } else {
    # Apply the residualize function to the entire dataset without considering waves.
    df <- residualise(data, age_residualise, sexres)
    if (rand) {
      df <- df %>% select(hhidpn, phenotype = std_resid)
    } else {
      df <- df %>% select(HHID, PN, phenotype = std_resid)
    }
  }

  # Save the processed data when requested, and always return it so the
  # same function can be reused when constructing composite phenotypes.
  if (save) {
    if (is.null(name) || !nzchar(name)) {
      stop("A non-empty name must be supplied when save = TRUE.")
    }
    if (!dir.exists(outputDataDir)) {
      dir.create(outputDataDir, recursive = TRUE)
    }
    fwrite(df, file.path(outputDataDir, paste0(name, "_", ancestry, ".pheno")))
  }

  return(df)
}


# Example call to the function
residualise.average.save(data=formerSmoker, average=FALSE, age_residualise=FALSE, name="SMCESS", rand=TRUE)

#########################################################
########## Construct phenotype: Alzheimer's ############
########################################################

# Transform data for residualising
ALZ <- data %>%
  # Selecting relevant columns and starting transformations
  select(hhidpn, contains("agem_m"), contains("rabyear"), contains("rabmonth"), contains("ragender"), contains("alzhe"), contains("alzhef")) %>%
  select(hhidpn, starts_with("r")) %>%
  rename(r10alzhef = raalzheef) %>%
  # Ensure all relevant columns are properly named and selected for transformation
  to_long(keys = "wave", values = c("alzhe", "alzhef", "age"),
          c(paste0("r", 10:14, "alzhe")),
          c(paste0("r", 10:14, "alzhef")),
          c(paste0("r", 10:14, "agem_m"))) %>%
  # Convert 'wave' to a clean numeric identifier
  mutate(wave = gsub("r", "", gsub("alzhe", "", wave)),
         alz = ifelse(alzhe == 1, 1, ifelse(age >= 65 & alzhe == 0, 0, NA)),
         alz = ifelse(alzhef == 0, alz, NA),
         # Convert gender to a binary numeric variable
         male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA))) %>%
  # Ensuring rabyear and rabmonth are numeric before creating 'dob'
  mutate(
    rabyear = as.numeric(rabyear),
    rabmonth = as.numeric(rabmonth),
    dob = ifelse(!is.na(rabyear) & !is.na(rabmonth), rabyear + (rabmonth / 12), NA_real_),
    dob2 = dob^2,
    male_dob = male * dob,
    male_dob2 = male * dob2
  ) %>%
  mutate(
    hhid = as.numeric(substr(hhidpn, 1, 5))  # Extract the first 5 digits and convert to numeric
  ) %>%
  # Filter out entries with essential missing values
  select(hhidpn, hhid, dob, dob2, male, male_dob, male_dob2, alz) %>%
  group_by(hhidpn) %>%
  # Aggregate data to get one record per hhidpn with pheno calculation
  mutate(pheno = ifelse(mean(alz, na.rm = TRUE) > 0, 1, ifelse(mean(alz, na.rm = TRUE) == 0, 0, NA)), 
         rn = row_number()) %>%
  ungroup() %>%
  # Filter to keep only the first record per hhidpn
  filter(rn == 1) %>%
  drop_na()
print(head(ALZ))  

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
ALZ <- ALZ %>%
  filter(hhid %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- ALZ %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(hhidpn)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in ALZ:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(ALZ)) 

# residualise, average, save
residualise.average.save(data=ALZ, average=F, age_residualise=F, rand = T, name="ALZ") 

########################################################
############### Construct phenotype: BMI ###############
########################################################

# Transform data for residualising
BMI <- data %>%
  select(hhidpn, contains("agem_m"), contains("gender"), contains("bmi")) %>%
  select(hhidpn, starts_with("r")) %>%
  to_long(keys = "wave", values = c("bmi", "age"),
                                  c(paste0("r", 1:14, "bmi")),
                                  c(paste0("r", 1:14, "agem_m"))) %>%
  mutate(
    age = as.numeric(age),  # Convert 'age' to numeric, handling possible NAs from conversion
    bmi = as.numeric(bmi),  # Convert 'bmi' to numeric, ensuring only valid numeric conversions are retained
    age2 = if_else(!is.na(age), age^2, NA_real_),  # Only calculate age2 if age is not NA
    male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA_real_)),  # Convert gender information with NA handling
    male_age = if_else(!is.na(male) & !is.na(age), male * age, NA_real_),
    male_age2 = if_else(!is.na(male_age), male_age * age, NA_real_),
    wave = as.numeric(gsub("r", "", gsub("bmi", "", wave)))  # Convert wave identifiers to numeric
  ) %>%
  select(-c(ragender, respagem_m)) %>%
  select(-contains("pmbmi")) %>%
  rename(pheno = bmi) %>%
  filter(!is.na(pheno) & !is.na(age) & !is.na(age2) & !is.na(male) & !is.na(male_age) & !is.na(male_age2))  %>% # Explicitly filter out NAs
  filter(pheno >= 10) # biologically implausible to have bmi <10

# Check for any non-numeric or infinite values before proceeding
BMI <- BMI %>%
  filter(!is.infinite(pheno) & !is.nan(pheno))  # Remove infinite and NaN values from 'pheno'

print("Transformed BMI data:")
print(head(BMI))


# Calculate summary statistics for the 'pheno' column
summary_stats <- BMI %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(hhidpn)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in BMI:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(BMI)) 

# residualise, average, save
residualise.average.save(data=BMI, average=T, age_residualise=T, rand = T, name="BMI")

########################################################
############### Construct phenotype: DPW ###############
########################################################

# transform data for residualising
DPW <- data %>%
  select(
    hhidpn,
    contains("agem_m"),
    contains("gender"),
    contains("drinkd"),
    contains("drinkn")
  ) %>%
  select(hhidpn, starts_with("r"), -"r1agem_m", -"r2agem_m") %>%
  to_long(keys = "wave", values = c("DrinksPerDay", "DaysDrinkingPerWeek", "age"),
                                  c(paste0("r", 3:14, "drinkn")),
                                  c(paste0("r", 3:14, "drinkd")),
                                  c(paste0("r", 3:14, "agem_m"))) %>%
  mutate(
    age2 = age^2,
    male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
    male_age = male * age,
    male_age2 = male * age2,
    wave = gsub("r", "", gsub("dpw", "", wave)),
    pheno = DrinksPerDay * DaysDrinkingPerWeek
  ) %>%
  select(-c(ragender, respagem_m, DrinksPerDay, DaysDrinkingPerWeek)) %>%
  drop_na()
  print(head(DPW))

# residualise, average, save
residualise.average.save(data=DPW, average=T, age_residualise=T, rand = T, name="DPW")


########################################################
## Construct phenotype: delaydisc (financial horizon) ##
########################################################

# transform data for residualising
delaydisc <- data %>%
  select(hhidpn,
         contains("agem_m"),
         contains("gender"),
         contains("finpln")) %>%
  select(hhidpn, starts_with("r")) %>%
  select(-c(r2agem_m, r3agem_m, r9agem_m, r10agem_m, respagem_m)) %>%
  to_long(keys = "wave", values = c("finpln", "age"),
                                  c(paste0("r", c(1, 4:8, 11:14), "finpln")),
                                  c(paste0("r", c(1, 4:8, 11:14), "agem_m"))) %>%
  mutate(
    age2 = age^2,
    male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
    male_age = male * age,
    male_age2 = male * age2,
    wave = gsub("r", "", gsub("finpln", "", wave)),
    pheno = case_when(
      finpln == "1.next few months" ~ 1,
      finpln == "2.next year" ~ 2,
      finpln == "3.next few years" ~ 3,
      finpln == "4.next 5-10 years" ~ 4,
      finpln == "5.longer than 10 years" ~ 5
    )
  ) %>%
  select(-c(finpln, ragender)) %>%
  drop_na()

# residualise, average, save
#residualise.average.save(data=delaydisc, average=T, age_residualise=T, rand = T, name="DELAYDISC")


########################################################
############## Construct phenotype: height #############
########################################################

# Transform data for residualising
height <- data %>%
  select(hhidpn,
         contains("agem_m"),
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("height")) %>%
  select(hhidpn, starts_with("r")) %>%
  to_long(keys = "wave", values = c("height", "age"),
                                  c(paste0("r", 1:14, "height")),
                                  c(paste0("r", 1:14, "agem_m"))) %>%
  mutate(
    height = as.numeric(height),  # Ensure height is numeric
    age = as.numeric(age),  # Ensure age is numeric
    rabyear = as.numeric(rabyear),  # Ensure rabyear is numeric
    rabmonth = as.numeric(rabmonth),  # Ensure rabmonth is numeric
    wave = gsub("r", "", gsub("height", "", wave)),
    male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
    dob = rabyear + (rabmonth / 12),
    dob2 = dob^2,
    male_dob = male * dob,
    male_dob2 = male * dob2
  ) %>%
  select(-c(ragender, respagem_m)) %>%
  group_by(hhidpn) %>%
  mutate(pheno = mean(height, na.rm=T), rn=row_number()) %>%
  ungroup() %>%
  filter(rn==1) %>%
  select(hhidpn, dob, dob2, male, male_dob, male_dob2, pheno) %>%
  drop_na()
  print(head(height))

# residualise, average, save
residualise.average.save(data=height, average=F, age_residualise=F, rand = T, name="HEIGHT")

##################################################################
############ Construct phenotype: smoking - eversmoke ############
##################################################################

# Transform data for residualising
smokerPhenotypes <- data %>%
  select(
    hhidpn,
    contains("rabyear"),
    contains("rabmonth"),
    contains("gender"),
    matches("r\\d+smokev")  # Selecting only the smokev columns
  ) %>%
  # Convert all smokev columns and date-related columns to numeric, handling non-numerics as NA
  mutate(
    across(matches("smokev"), ~ifelse(is.na(as.numeric(.)), NA_real_, as.numeric(.))),
    rabyear = as.numeric(rabyear),
    rabmonth = as.numeric(rabmonth),
    male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
    dob = ifelse(is.na(rabyear) | is.na(rabmonth), NA_real_, rabyear + (rabmonth / 12)),  # Safely calculate dob
    dob2 = dob^2,
    male_dob = male * dob,
    male_dob2 = male * dob2
  ) %>%
  rowwise() %>%
  mutate(
    # Calculate everSmoker based on smokev columns
    pheno = if (all(is.na(c_across(matches("smokev"))))) {
      NA_real_  # Set to NA if all values are missing
    } else if (any(c_across(matches("smokev")) == 1, na.rm = TRUE)) {
      1  # Set to 1 if any wave reports smoking
    } else {
      0  # Set to 0 if all non-missing values are 0
    }
  ) %>%
  mutate(
    hhid = as.numeric(substr(hhidpn, 1, 5))  # Extract the first 5 digits and convert to numeric
  ) %>%
  select(hhidpn, hhid, dob, dob2, male, male_dob, male_dob2, pheno) %>%
  ungroup() %>%
  drop_na()
print(head(smokerPhenotypes, 10))  

# Merge the new phenotype back with the original data
diagnostic_data <- data %>%
  left_join(smokerPhenotypes, by = "hhidpn") %>%
  select(hhidpn, matches("r\\d+smokev"), dob, dob2, male, male_dob, male_dob2, pheno)

# Print the structure and a sample of the merged data for diagnostic purposes
print("Structure of Diagnostic Data with everSmoker:")
str(diagnostic_data)

print("First 10 Rows for Diagnostic Check with everSmoker:")
print(head(diagnostic_data, 10))

print("Last 10 Rows for Diagnostic Check with everSmoker:")
print(tail(diagnostic_data, 10))

# Filter the smokerPhenotypes data frame for individuals with HHID in the requested ancestry HHIDs list
smokerPhenotypes <- smokerPhenotypes %>%
  filter(hhid %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- smokerPhenotypes %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(hhidpn)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in smokerPhenotypes:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(smokerPhenotypes)) 

# residualise, average, save
residualise.average.save(data=smokerPhenotypes, average=F, age_residualise=F, rand = T, name="EVERSMOKE")

#########################################################################################
############ Construct phenotype: smoking - former smoker/ smoking cessation ############
#########################################################################################

# transform data for residualising
formerSmoker <- data %>%
  select(hhidpn,
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("smoken"),
         contains("smokev")) %>%
  select(hhidpn, starts_with("r")) %>%
  to_long(keys = "wave", values = c("smokev","smoken"),
                                  c(paste0("r", 1:14, "smokev")),
                                  c(paste0("r", 1:14, "smoken"))) %>%
  mutate(smokev = ifelse(smokev == 1, 1, ifelse(smokev == 0, 0, NA)),
         smoken = ifelse(smoken == 1, 1, ifelse(smoken == 0, 0, NA)),
         smokecess = ifelse(smokev == 1 & smoken == 0, 1, ifelse(smoken == 1, 0, NA))) %>%
  group_by(hhidpn) %>%
  mutate(pheno = ifelse(mean(smokecess, na.rm = TRUE) > 0, 1, ifelse(mean(smokecess, na.rm = TRUE) == 0, 0, NA)),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn == 1) %>%
  mutate(rabyear = as.numeric(rabyear),  # Convert rabyear to numeric
         rabmonth = as.numeric(rabmonth),  # Convert rabmonth to numeric
         male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
         dob = rabyear + (rabmonth / 12),  # Calculate dob
         dob2 = dob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2) %>%
  mutate(
    hhid = as.numeric(substr(hhidpn, 1, 5))  # Extract the first 5 digits and convert to numeric
  ) %>%
  select(hhidpn, hhid, dob, dob2, male, male_dob, male_dob2, pheno) %>%
  drop_na()

# Filter the formerSmoker data frame for individuals with HHID in the requested ancestry HHIDs list
formerSmoker <- formerSmoker %>%
  filter(hhid %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- formerSmoker %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(hhidpn)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in formerSmoker:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(formerSmoker)) 

# residualise, average, save
residualise.average.save(data=formerSmoker, average=F, age_residualise=F, rand = T, name="SMCESS")

#################################################################################
################ Construct phenotype: Educational Attainement (EA) ###############
##################################################################################

# Transform data for residualising
EA <- data %>%
  select(hhidpn,
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("raedyrs")) %>%
  select(hhidpn, starts_with("r")) %>%
  mutate(
    across(c(rabyear, rabmonth, raedyrs), ~ as.numeric(.)),  # Ensure columns are numeric
    male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
    dob = rabyear + (rabmonth / 12),
    dob2 = dob^2,
    male_dob = male * dob,
    male_dob2 = male * dob2,
  ) %>%
  mutate(
    hhid = as.numeric(substr(hhidpn, 1, 5))  # Extract the first 5 digits and convert to numeric
  ) %>%
  select(hhidpn, hhid, dob, dob2, male, male_dob, male_dob2, pheno = raedyrs) %>%
  drop_na()
print(head(EA$hhid))

# Filter EA based on HHIDs
EA <- EA %>%
  filter(hhid %in% ancestry_hhids)

# Calculate summary statistics for the 'pheno' column
summary_stats <- EA %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(hhidpn)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in EA:")
print(summary_stats)

# Print the number of records after filtering for requested ancestry
print(sprintf("Number of records after filtering: %d", nrow(EA)))

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(EA)) 

# residualise, average, save
residualise.average.save(data=EA, average=F, age_residualise=F, rand = T, name="EA")

########################################################
############ Construct phenotype: loneliness ###########
########################################################

# Convert haven_labelled columns to numeric
data <- data %>%
  purrr::map_df(~ if(is.labelled(.x)) as.numeric(.x) else .x)

# Ensure all flone columns are of the same type
data <- data %>%
  mutate(across(starts_with("r") & contains("flone"), as.character))

# Transform data for residualising
lonely <- data %>%
  select(hhidpn, contains("agem_m"), contains("gender"), contains("flone")) %>%
  mutate(r1flone = r1flonex) %>% # Add the missing variable to avoid errors in select
  select(hhidpn,
         starts_with("r"),
         r1flone) %>% # Select the renamed variable

  # First pivot for "flone"
  pivot_longer(cols = starts_with("r") & contains("flone"),
               names_to = "wave",
               values_to = "flone") %>%

  # Ensure the wave columns are of the same type
  mutate(across(starts_with("r") & contains("agem_m"), as.numeric)) %>%
  
  # Second pivot for "age"
  pivot_longer(cols = starts_with("r") & contains("agem_m"),
               names_to = "wave_age",
               values_to = "age") %>%

  mutate(wave = gsub("flone", "", wave),
         age2 = age^2,
         male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
         male_age = male * age,
         male_age2 = male * age2,
         pheno = ifelse(flone == 1, 1, ifelse(flone == 2, .67, ifelse(flone == 3, .33, ifelse(flone == 0, 0, NA)))),
         wave = gsub("r", "", gsub("flone", "", wave))) %>%
  select(-c(ragender)) %>%
  select(-flone) %>%
  drop_na()

# Print the first few rows to check the transformation
head(lonely)

# residualise, average, save
residualise.average.save(data=lonely, average=T, age_residualise=T, rand = T, name="LONELY")


########################################################
########## Construct phenotype: NEBmen, NEBwom #########
########################################################

# Convert haven_labelled columns to numeric
data <- data %>%
  purrr::map_df(~ if(is.labelled(.x)) as.numeric(.x) else .x)

# Ensure all relevant columns are of the appropriate type
data <- data %>%
  mutate(across(starts_with("r") & contains("agem_m"), as.character),
         across(c(rabyear, rabmonth), as.numeric))

# Transform data for residualising
neb <- data %>%
  select(hhidpn, contains("raevbrn"), contains("gender"), contains("agem_m")) %>%
  mutate(male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA))) %>%
  
  # Pivot for "age"
  pivot_longer(cols = starts_with("r") & contains("agem_m"),
               names_to = "wave_age",
               values_to = "age_m") %>%
  mutate(wave = sub("agem_m", "", wave_age),
         age = as.numeric(age_m) / 12) %>%
  group_by(hhidpn, male) %>%
  summarize(min_age = if(all(is.na(age))) NA else min(age, na.rm = TRUE)) %>%
  filter((male == 0 & min_age >= 45) | (male == 1 & min_age >= 55)) %>%
  ungroup() %>%
  
  # Join back to original data to get other variables
  left_join(data, by = "hhidpn") %>%
  mutate(dob = rabyear + (rabmonth / 12),
         dob2 = dob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = as.numeric(raevbrn)) %>%  # Convert pheno to numeric
  mutate(
    hhid = as.numeric(substr(hhidpn, 1, 5))  # Extract the first 5 digits and convert to numeric
  ) %>%
  select(hhidpn, hhid, dob, dob2, male, male_dob, male_dob2, min_age, pheno) %>%
  drop_na() %>%
  filter_all(all_vars(!is.na(.) & !is.nan(.) & !is.infinite(.)))  # Ensure no NA/NaN/Inf values

# Filter out non-numeric values in pheno
neb <- neb %>%
  filter(!is.na(pheno))

# Diagnostic for age-based filtering 
min_age_female <- min(neb$min_age[neb$male == 0], na.rm = TRUE)
print(paste("Minimum age for females:", min_age_female))
min_age_male <- min(neb$min_age[neb$male == 1], na.rm = TRUE)
print(paste("Minimum age for males:", min_age_male))
# Count of males and females after filtering
print(table(neb$male))
# Summary of minimum age for females
print(summary(neb$min_age[neb$male == 0]))
# First few rows of female records after filtering
print(head(neb[neb$male == 0, ]))

# Restrict by gender to get NEBmen and NEBwomen stats 
neb <- neb %>%
  filter(hhid %in% ancestry_hhids, male == 1)

# Calculate summary statistics for the 'pheno' column
summary_stats <- neb %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(hhidpn)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in NEB:")
print(summary_stats)

# Residualise, average, save
residualise.average.save(data=neb, average=F, age_residualise=F, name="NEB", rand = T, neb=T)

########################################################
############### Construct phenotype: risk ##############
########################################################

# transform data for residualising
risk <- data %>%
  select(hhidpn, contains("agem_m"), contains("gender"), contains("risk")) %>%
  select(hhidpn,
         starts_with("r"),
         -c(r2agem_m, r3agem_m, r9agem_m, r10agem_m, r11agem_m, r12agem_m, r13agem_m)) %>%
  gather(key="wave", value="value", paste0("r", c(1, 4:8), "risk"), paste0("r", 4:8, "risk6"), paste0("r", c(1, 4:8), "agem_m")) %>%
  mutate(wave = recode(wave,

                       "r1agem_m"="r1_agem",
                       "r4agem_m"="r4_agem",
                       "r5agem_m"="r5_agem",
                       "r6agem_m"="r6_agem",
                       "r7agem_m"="r7_agem",
                       "r8agem_m"="r8_agem",

                       "r1risk" = "r1_risk",
                       "r4risk" = "r4_risk",
                       "r5risk" = "r5_risk",
                       "r6risk" = "r6_risk",
                       "r7risk" = "r7_risk",
                       "r8risk" = "r8_risk",

                       "r4risk6" = "r4_risk6",
                       "r5risk6" = "r5_risk6",
                       "r6risk6" = "r6_risk6",
                       "r7risk6" = "r7_risk6",
                       "r8risk6" = "r8_risk6")) %>%
  separate("wave", c("wave", "var")) %>%
  mutate(wave = substr(wave, 2, 2)) %>%
  spread("var", "value") %>%
  mutate(age = as.numeric(agem),
         age2 = age^2,
         male = ifelse(ragender == 1, 1, ifelse(ragender == 2, 0, NA)),
         male_age = male * age,
         male_age2 = male * age2,
         risk = as.numeric(substr(risk, 1, 1)),
         risk6 = as.numeric(substr(risk6, 1, 1)),
         pheno = case_when(!is.na(risk6) ~ risk6 * 4 / 6,
                           !is.na(risk)  ~ risk)) %>%
  drop_na(age, male, pheno)

# residualise, average, save
residualise.average.save(data=risk, average=T, age_residualise=T, name="RISK", rand = T, neb=F)


########################################################
################## List HRS data paths #################
########################################################

path_1992 <- file.path(inputDataDir, "1992_merged_respondent_level.txt")
path_1993 <- file.path(inputDataDir, "1993_merged_respondent_level.txt")
path_1994 <- file.path(inputDataDir, "1994_merged_respondent_level.txt")
path_1995 <- file.path(inputDataDir, "1995_merged_respondent_level.txt")
path_1996 <- file.path(inputDataDir, "1996_merged_respondent_level.txt")
path_1998 <- file.path(inputDataDir, "1998_merged_respondent_level.txt")
path_2000 <- file.path(inputDataDir, "2000_merged_respondent_level.txt")
path_2002 <- file.path(inputDataDir, "2002_merged_respondent_level.txt")
path_2004 <- file.path(inputDataDir, "HRS_2004_data.txt")
path_2006 <- file.path(inputDataDir, "HRS_2006_data.txt")
path_2008 <- file.path(inputDataDir, "HRS_2008_data.txt")
path_2010 <- file.path(inputDataDir, "HRS_2010_data.txt")
path_2012 <- file.path(inputDataDir, "HRS_2012_data.txt")
path_2014 <- file.path(inputDataDir, "HRS_2014_data.txt")
path_2016 <- file.path(inputDataDir, "HRS_2016_data.txt")
path_2018 <- file.path(inputDataDir, "HRS_2018_data.txt")
path_tracker <- file.path(inputDataDir, "TRK2018TR_R_stataversion12.dta")


########################################################
############ Load and wrangle HRS wave data ############
########################################################

data_tracker <- read.dta(path_tracker) %>%
  mutate(HHID = as.numeric(HHID),
         PN = as.numeric(PN),
         sex = as.numeric(GENDER),
         mob = as.numeric(BIRTHMO),
         yob = as.numeric(BIRTHYR)) %>%
  select(HHID, PN, sex, mob, yob) %>%
  filter(mob != 0 & yob != 0)

data_1992 <- fread(path_1992) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         age_1992=V46,
         b_moi_1992=V23, b_yoi_1992=V25, e_moi_1992=V26, e_yoi_1992=V28,
         moi_1992 = e_moi_1992, yoi_1992 = e_yoi_1992)

data_1993 <- fread(path_1993) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         age_1993=AGE,
         moi_1993=V359, yoi_1993=V361)

data_1994 <- fread(path_1994) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         age_1994=W104,
         f_moi_1994=W120, f_yoi_1994=W122, b_moi_1994=W56, b_yoi_1994=W58, e_moi_1994=W59, e_yoi_1994=W61,
         moi_1994 = e_moi_1994, yoi_1994 = e_yoi_1994)

data_1995 <- fread(path_1995) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         c_moi_1995=D391, f_moi_1995=D398, moi_1995 = c_moi_1995,
         yoi_1995=D393)

data_1996 <- fread(path_1996) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_1996=E391, yoi_1996=E393)

data_1998 <- fread(path_1998) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_1998=F704, yoi_1998=F703,
         e_moi_1998=F697, e_yoi_1998=F699)

data_2000 <- fread(path_2000) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2000=G775, yoi_2000=G774,
         e_moi_2000=G768, e_yoi_2000=G770)

data_2002 <- fread(path_2002) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2002=HA500, yoi_2002=HA501)

data_2004 <- fread(path_2004) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2004=JA500, yoi_2004=JA501)

data_2006 <- fread(path_2006) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2006=KA500, yoi_2006=KA501)

data_2008 <- fread(path_2008) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2008=LA500, yoi_2008=LA501)

data_2010 <- fread(path_2010) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2010=MA500, yoi_2010=MA501)

data_2012 <- fread(path_2012) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2012=NA500, yoi_2012=NA501)

data_2014 <- fread(path_2014) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2014=OA500, yoi_2014=OA501)

data_2016 <- fread(path_2016) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2016=PA500, yoi_2016=PA501)

data_2018 <- fread(path_2018) %>%
  full_join(data_tracker, by = c("HHID", "PN")) %>%
  mutate(sex, mob, yob,
         moi_2018=QA500, yoi_2018=QA501)

########################################################
################ Construct phenotype: SWB ##############
########################################################

# Subjective well-being (SWB) combines positive affect (PA) and
# life satisfaction (LS)

force_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

reverse_swb_pa <- function(x) {
  dplyr::recode(x, `1` = 5, `2` = 4, `4` = 2, `5` = 1, .default = x)
}

row_mean_with_missing_limit <- function(data, variables, max_missing) {
  item_data <- data %>%
    select(all_of(variables))

  missing_count <- rowSums(is.na(item_data))
  score <- rowMeans(item_data, na.rm = TRUE)
  score[is.nan(score) | missing_count >= max_missing] <- NA_real_
  score
}

construct_swb_wave <- function(data_wave, year, sex_var, mob_var, yob_var,
                               moi_var, yoi_var, pa_vars, ls_vars,
                               pa_max_missing, ls_max_missing = 3) {

  required_vars <- c(
    "HHID", "PN", sex_var, mob_var, yob_var, moi_var, yoi_var,
    pa_vars, ls_vars
  )

  missing_vars <- setdiff(required_vars, names(data_wave))
  if (length(missing_vars) > 0) {
    stop(
      "Missing variables for SWB ", year, ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  wave_data <- data_wave %>%
    select(all_of(required_vars)) %>%
    mutate(
      across(all_of(c(sex_var, mob_var, yob_var, moi_var, yoi_var,
                      pa_vars, ls_vars)), force_numeric),
      across(all_of(pa_vars), reverse_swb_pa)
    )

  wave_data$PA <- row_mean_with_missing_limit(
    wave_data, pa_vars, max_missing = pa_max_missing
  )
  wave_data$LS <- row_mean_with_missing_limit(
    wave_data, ls_vars, max_missing = ls_max_missing
  )

  wave_data %>%
    transmute(
      HHID,
      PN,
      wave = year,
      PA,
      LS,
      sex = .data[[sex_var]],
      mob = .data[[mob_var]],
      yob = .data[[yob_var]],
      moi = .data[[moi_var]],
      yoi = .data[[yoi_var]]
    ) %>%
    filter(mob %in% 1:12, moi %in% 1:12) %>%
    mutate(
      dob = yob + mob / 12,
      doi = yoi + moi / 12,
      age = doi - dob,
      age2 = age^2,
      male = 2 - sex,
      male_age = male * age,
      male_age2 = male * age2
    )
}

# Construct PA and LS scales separately in each wave.
swb_waves <- bind_rows(
  construct_swb_wave(
    data_2006, 2006,
    sex_var = "KX060_R", mob_var = "KX004_R", yob_var = "KX067_R",
    moi_var = "KA500", yoi_var = "KA501",
    pa_vars = paste0("KLB027", LETTERS[1:8]),
    ls_vars = paste0("KLB003", LETTERS[1:5]),
    pa_max_missing = 5
  ),
  construct_swb_wave(
    data_2008, 2008,
    sex_var = "LX060_R", mob_var = "LX004_R", yob_var = "LX067_R",
    moi_var = "LA500", yoi_var = "LA501",
    pa_vars = paste0("LLB027", c("C", "D", "F", "G", "H", "K", "P", "Q", "T", "U", "V", "X", "Y")),
    ls_vars = paste0("LLB003", LETTERS[1:5]),
    pa_max_missing = 7
  ),
  construct_swb_wave(
    data_2010, 2010,
    sex_var = "MX060_R", mob_var = "MX004_R", yob_var = "MX067_R",
    moi_var = "MA500", yoi_var = "MA501",
    pa_vars = paste0("MLB027", c("C", "D", "F", "G", "H", "K", "P", "Q", "T", "U", "V", "X", "Y")),
    ls_vars = paste0("MLB003", LETTERS[1:5]),
    pa_max_missing = 7
  ),
  construct_swb_wave(
    data_2012, 2012,
    sex_var = "NX060_R", mob_var = "NX004_R", yob_var = "NX067_R",
    moi_var = "NA500", yoi_var = "NA501",
    pa_vars = paste0("NLB027", c("C", "D", "F", "G", "H", "K", "P", "Q", "T", "U", "V", "X", "Y")),
    ls_vars = paste0("NLB003", LETTERS[1:5]),
    pa_max_missing = 7
  ),
  construct_swb_wave(
    data_2014, 2014,
    sex_var = "OX060_R", mob_var = "OX004_R", yob_var = "OX067_R",
    moi_var = "OA500", yoi_var = "OA501",
    pa_vars = paste0("OLB026", c("C", "D", "F", "G", "H", "K", "P", "Q", "T", "U", "V", "X", "Y")),
    ls_vars = paste0("OLB002", LETTERS[1:5]),
    pa_max_missing = 7
  ),
  construct_swb_wave(
    data_2016, 2016,
    sex_var = "PX060_R", mob_var = "PX004_R", yob_var = "PX067_R",
    moi_var = "PA500", yoi_var = "PA501",
    pa_vars = paste0("PLB026", c("C", "D", "F", "G", "H", "K", "P", "Q", "T", "U", "V", "X", "Y")),
    ls_vars = paste0("PLB002", LETTERS[1:5]),
    pa_max_missing = 7
  )
)

PA <- swb_waves %>%
  filter(!is.na(PA)) %>%
  transmute(HHID, PN, wave, pheno = PA,
            age, age2, male, male_age, male_age2) %>%
  drop_na()

LS <- swb_waves %>%
  filter(!is.na(LS)) %>%
  transmute(HHID, PN, wave, pheno = LS,
            age, age2, male, male_age, male_age2) %>%
  drop_na()

# Residualise separately by wave and average each component across waves
# using the same shared function as the other longitudinal phenotypes.
SWBPA <- residualise.average.save(
  data = PA,
  average = TRUE,
  age_residualise = TRUE,
  sexres = TRUE,
  rand = FALSE,
  save = FALSE
) %>%
  rename(PA = phenotype)

SWBLS <- residualise.average.save(
  data = LS,
  average = TRUE,
  age_residualise = TRUE,
  sexres = TRUE,
  rand = FALSE,
  save = FALSE
) %>%
  rename(LS = phenotype)

# Preserve the original implementation: SWB is the sum, not the arithmetic mean,
# of the person-level PA and LS scores, and requires both components.
SWB <- inner_join(SWBPA, SWBLS, by = c("HHID", "PN")) %>%
  filter(HHID %in% ancestry_hhids) %>%
  mutate(phenotype = PA + LS) %>%
  select(HHID, PN, phenotype)

fwrite(SWB, file.path(outputDataDir, paste0("SWB_", ancestry, ".pheno")))


########################################################
############### Construct phenotype: CIDI ##############
########################################################

# Construct the CIDI Short Form depression phenotype directly in R.
# This follows 10.1.5_construct_HRS_CIDI.do without temporary Stata files,
# while including the available 2014 and 2016 scores in the probability average.

cidi_years <- c(1995, 1996, seq(1998, 2016, 2))
cidi_probability_years <- cidi_years

# Raw item mappings used in each HRS wave.
cidi_item_map <- list(
  `1995` = c(
    DEPR_screen = "D1006", ANH_screen = "D1028",
    DEPR_intensity = "D1007", DEPR_freq = "D1008",
    ANH_intensity = "D1029", ANH_freq = "D1030",
    DEPR_anh = "D1009", DEPR_tired = "D1010",
    DEPR_lossApp = "D1011", DEPR_incrApp = "D1012",
    DEPR_sleep = "D1013", DEPR_sleepFreq = "D1014",
    DEPR_ctrate = "D1015", DEPR_down = "D1016",
    DEPR_death = "D1017",
    ANH_tired = "D1031", ANH_lossApp = "D1032",
    ANH_incrApp = "D1033", ANH_sleep = "D1034",
    ANH_sleepFreq = "D1035", ANH_ctrate = "D1036",
    ANH_down = "D1037", ANH_death = "D1038"
  ),
  `1996` = c(
    DEPR_screen = "E1006", ANH_screen = "E1028",
    DEPR_intensity = "E1007", DEPR_freq = "E1008",
    ANH_intensity = "E1029", ANH_freq = "E1030",
    DEPR_anh = "E1009", DEPR_tired = "E1010",
    DEPR_lossApp = "E1011", DEPR_incrApp = "E1012",
    DEPR_sleep = "E1013", DEPR_sleepFreq = "E1014",
    DEPR_ctrate = "E1015", DEPR_down = "E1016",
    DEPR_death = "E1018",
    ANH_tired = "E1031", ANH_lossApp = "E1032",
    ANH_incrApp = "E1033", ANH_sleep = "E1034",
    ANH_sleepFreq = "E1035", ANH_ctrate = "E1036",
    ANH_down = "E1037", ANH_death = "E1038"
  ),
  `1998` = c(
    DEPR_screen = "F1323", ANH_screen = "F1345",
    DEPR_intensity = "F1324", DEPR_freq = "F1325",
    ANH_intensity = "F1346", ANH_freq = "F1347",
    DEPR_anh = "F1326", DEPR_tired = "F1327",
    DEPR_lossApp = "F1328", DEPR_incrApp = "F1329",
    DEPR_sleep = "F1330", DEPR_sleepFreq = "F1331",
    DEPR_ctrate = "F1332", DEPR_down = "F1333",
    DEPR_death = "F1334",
    ANH_tired = "F1348", ANH_lossApp = "F1349",
    ANH_incrApp = "F1350", ANH_sleep = "F1351",
    ANH_sleepFreq = "F1352", ANH_ctrate = "F1353",
    ANH_down = "F1354", ANH_death = "F1355"
  ),
  `2000` = c(
    DEPR_screen = "G1456", ANH_screen = "G1478",
    DEPR_intensity = "G1457", DEPR_freq = "G1458",
    ANH_intensity = "G1479", ANH_freq = "G1480",
    DEPR_anh = "G1459", DEPR_tired = "G1460",
    DEPR_lossApp = "G1461", DEPR_incrApp = "G1462",
    DEPR_sleep = "G1463", DEPR_sleepFreq = "G1464",
    DEPR_ctrate = "G1465", DEPR_down = "G1466",
    DEPR_death = "G1467",
    ANH_tired = "G1481", ANH_lossApp = "G1482",
    ANH_incrApp = "G1483", ANH_sleep = "G1484",
    ANH_sleepFreq = "G1485", ANH_ctrate = "G1486",
    ANH_down = "G1487", ANH_death = "G1488"
  )
)

# From 2002 onward the item positions are constant and only the wave prefix changes.
cidi_later_prefixes <- c(
  `2002` = "H", `2004` = "J", `2006` = "K", `2008` = "L",
  `2010` = "M", `2012` = "N", `2014` = "O", `2016` = "P"
)

for (year in names(cidi_later_prefixes)) {
  prefix <- cidi_later_prefixes[[year]]
  cidi_item_map[[year]] <- c(
    DEPR_screen = paste0(prefix, "C150"),
    DEPR_intensity = paste0(prefix, "C151"),
    DEPR_freq = paste0(prefix, "C152"),
    DEPR_anh = paste0(prefix, "C153"),
    DEPR_tired = paste0(prefix, "C154"),
    DEPR_lossApp = paste0(prefix, "C155"),
    DEPR_incrApp = paste0(prefix, "C156"),
    DEPR_sleep = paste0(prefix, "C157"),
    DEPR_sleepFreq = paste0(prefix, "C158"),
    DEPR_ctrate = paste0(prefix, "C159"),
    DEPR_down = paste0(prefix, "C160"),
    DEPR_death = paste0(prefix, "C161"),
    ANH_screen = paste0(prefix, "C167"),
    ANH_intensity = paste0(prefix, "C168"),
    ANH_freq = paste0(prefix, "C169"),
    ANH_tired = paste0(prefix, "C170"),
    ANH_lossApp = paste0(prefix, "C171"),
    ANH_incrApp = paste0(prefix, "C172"),
    ANH_sleep = paste0(prefix, "C173"),
    ANH_sleepFreq = paste0(prefix, "C174"),
    ANH_ctrate = paste0(prefix, "C175"),
    ANH_down = paste0(prefix, "C176"),
    ANH_death = paste0(prefix, "C177")
  )
}

cidi_wave_data <- list(
  `1995` = data_1995, `1996` = data_1996, `1998` = data_1998,
  `2000` = data_2000, `2002` = data_2002, `2004` = data_2004,
  `2006` = data_2006, `2008` = data_2008, `2010` = data_2010,
  `2012` = data_2012, `2014` = data_2014, `2016` = data_2016
)

construct_cidi_wave <- function(data_wave, year, item_map) {
  condition_true <- function(x) !is.na(x) & x

  interview_vars <- c(
    `1995` = "moi_1995", `1996` = "moi_1996", `1998` = "moi_1998",
    `2000` = "moi_2000", `2002` = "moi_2002", `2004` = "moi_2004",
    `2006` = "moi_2006", `2008` = "moi_2008", `2010` = "moi_2010",
    `2012` = "moi_2012", `2014` = "moi_2014", `2016` = "moi_2016"
  )
  interview_year_vars <- c(
    `1995` = "yoi_1995", `1996` = "yoi_1996", `1998` = "yoi_1998",
    `2000` = "yoi_2000", `2002` = "yoi_2002", `2004` = "yoi_2004",
    `2006` = "yoi_2006", `2008` = "yoi_2008", `2010` = "yoi_2010",
    `2012` = "yoi_2012", `2014` = "yoi_2014", `2016` = "yoi_2016"
  )

  required_vars <- c(
    "HHID", "PN", "sex", "mob", "yob",
    interview_vars[[as.character(year)]],
    interview_year_vars[[as.character(year)]],
    unname(item_map)
  )
  missing_vars <- setdiff(required_vars, names(data_wave))

  if (length(missing_vars) > 0) {
    stop(
      "Missing variables for CIDI ", year, ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  moi_name <- interview_vars[[as.character(year)]]
  yoi_name <- interview_year_vars[[as.character(year)]]

  items <- data_wave %>%
    select(all_of(required_vars)) %>%
    rename(!!!setNames(unname(item_map), names(item_map))) %>%
    mutate(
      across(everything(), force_numeric),
      across(all_of(names(item_map)), ~ replace(.x, .x >= 7, NA_real_)),
      dob = if_else(
        mob %in% 1:12 & is.finite(yob),
        yob + mob / 12,
        NA_real_
      ),
      doi = if_else(
        .data[[moi_name]] %in% 1:12 & is.finite(.data[[yoi_name]]),
        .data[[yoi_name]] + .data[[moi_name]] / 12,
        NA_real_
      ),
      age = doi - dob,
      age = if_else(age >= 0 & age <= 120, age, NA_real_),
      male = case_when(
        sex == 1 ~ 1,
        sex == 2 ~ 0,
        TRUE ~ NA_real_
      )
    )

  # The anhedonia screen itself counts as one anhedonia symptom.
  symptoms_DEPR <- rep(0, nrow(items))
  symptoms_ANH <- rep(1, nrow(items))

  for (screen in c("DEPR", "ANH")) {
    symptom_count <- if (screen == "DEPR") symptoms_DEPR else symptoms_ANH

    symptom_count <- symptom_count +
      ifelse(condition_true(items[[paste0(screen, "_tired")]] == 1), 1, 0) +
      ifelse(
        condition_true(
          items[[paste0(screen, "_lossApp")]] == 1 |
            items[[paste0(screen, "_incrApp")]] == 1
        ),
        1, 0
      ) +
      ifelse(
        condition_true(
          items[[paste0(screen, "_sleep")]] == 1 &
            items[[paste0(screen, "_sleepFreq")]] <= 2
        ),
        1, 0
      ) +
      ifelse(condition_true(items[[paste0(screen, "_ctrate")]] == 1), 1, 0) +
      ifelse(condition_true(items[[paste0(screen, "_down")]] == 1), 1, 0) +
      ifelse(condition_true(items[[paste0(screen, "_death")]] == 1), 1, 0)

    if (screen == "DEPR") {
      symptoms_DEPR <- symptom_count
    } else {
      symptoms_ANH <- symptom_count
    }
  }

  symptoms_DEPR <- symptoms_DEPR + ifelse(condition_true(items$DEPR_anh == 1), 1, 0)

  cidi <- rep(NA_real_, nrow(items))

  # Respondents who do not pass either screen receive -1.
  cidi[
    condition_true(
    items$DEPR_screen == 5 & items$ANH_screen == 5
    )
  ] <- -1
  cidi[
    condition_true(
    items$DEPR_screen == 5 & items$ANH_screen == 1 &
      (items$ANH_intensity > 2 | items$ANH_freq > 2)
    )
  ] <- -1
  cidi[
    condition_true(
    items$DEPR_screen == 1 &
      (items$DEPR_intensity > 2 | items$DEPR_freq > 2) &
      items$ANH_screen == 5
    )
  ] <- -1
  cidi[
    condition_true(
    items$DEPR_screen == 1 &
      (items$DEPR_intensity > 2 | items$DEPR_freq > 2) &
      items$ANH_screen == 1 &
      (items$ANH_intensity > 2 | items$ANH_freq > 2)
    )
  ] <- -1

  # Assign the symptom count for respondents who pass a screen.
  cidi[
    condition_true(
    items$DEPR_screen == 1 &
      items$DEPR_intensity <= 2 & items$DEPR_freq <= 2
    )
  ] <- symptoms_DEPR[
    condition_true(
    items$DEPR_screen == 1 &
      items$DEPR_intensity <= 2 & items$DEPR_freq <= 2
    )
  ]

  cidi[
    condition_true(
    items$DEPR_screen == 5 & items$ANH_screen == 1 &
      items$ANH_intensity <= 2 & items$ANH_freq <= 2
    )
  ] <- symptoms_ANH[
    condition_true(
    items$DEPR_screen == 5 & items$ANH_screen == 1 &
      items$ANH_intensity <= 2 & items$ANH_freq <= 2
    )
  ]

  cidi[
    condition_true(
    items$DEPR_screen == 1 &
      (items$DEPR_intensity > 2 | items$DEPR_freq > 2) &
      items$ANH_screen == 1 &
      items$ANH_intensity <= 2 & items$ANH_freq <= 2
    )
  ] <- symptoms_ANH[
    condition_true(
    items$DEPR_screen == 1 &
      (items$DEPR_intensity > 2 | items$DEPR_freq > 2) &
      items$ANH_screen == 1 &
      items$ANH_intensity <= 2 & items$ANH_freq <= 2
    )
  ]

  # A passing screen must have complete information for all required symptoms.
  depr_pass <- condition_true(
    items$DEPR_screen == 1 &
      items$DEPR_intensity <= 2 & items$DEPR_freq <= 2
  )
  depr_incomplete <-
    is.na(items$DEPR_anh) |
    is.na(items$DEPR_tired) |
    (is.na(items$DEPR_lossApp) & is.na(items$DEPR_incrApp)) |
    is.na(items$DEPR_sleep) |
    condition_true(items$DEPR_sleep == 1) & is.na(items$DEPR_sleepFreq) |
    is.na(items$DEPR_ctrate) |
    is.na(items$DEPR_down) |
    is.na(items$DEPR_death)

  anh_pass <- condition_true(
    items$ANH_screen == 1 &
      items$ANH_intensity <= 2 & items$ANH_freq <= 2
  )
  anh_incomplete <-
    is.na(items$ANH_tired) |
    (is.na(items$ANH_lossApp) & is.na(items$ANH_incrApp)) |
    is.na(items$ANH_sleep) |
    condition_true(items$ANH_sleep == 1) & is.na(items$ANH_sleepFreq) |
    is.na(items$ANH_ctrate) |
    is.na(items$ANH_down) |
    is.na(items$ANH_death)

  cidi[depr_pass & depr_incomplete] <- NA_real_
  cidi[anh_pass & anh_incomplete] <- NA_real_

  tibble(
    HHID = items$HHID,
    PN = items$PN,
    wave = as.integer(year),
    CIDI = cidi,
    age = items$age,
    male = items$male
  )
}

CIDI_long <- purrr::imap_dfr(
  cidi_wave_data,
  ~ construct_cidi_wave(.x, .y, cidi_item_map[[.y]])
)

CIDI_wide <- CIDI_long %>%
  select(HHID, PN, wave, CIDI, age, male) %>%
  group_by(HHID, PN) %>%
  mutate(
    # Sex comes from the tracker and should be constant across waves.
    male = dplyr::first(male[!is.na(male)], default = NA_real_)
  ) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(HHID, PN, male),
    names_from = wave,
    values_from = c(CIDI, age),
    names_glue = "{.value}_{wave}"
  ) %>%
  mutate(
    CIDIcount = rowSums(
      !is.na(across(all_of(paste0("CIDI_", cidi_years))))
    )
  ) %>%
  filter(CIDIcount > 0)

cidi_probability <- c(
  `-1` = 0,
  `0` = 0.0001,
  `1` = 0.0568,
  `2` = 0.2352,
  `3` = 0.5542,
  `4` = 0.8125,
  `5` = 0.8895,
  `6` = 0.8895,
  `7` = 0.9083
)

for (year in cidi_probability_years) {
  cidi_var <- paste0("CIDI_", year)
  prob_var <- paste0("Depr_prob_", year)
  age_var <- paste0("age_", year)

  CIDI_wide[[prob_var]] <- unname(
    cidi_probability[as.character(CIDI_wide[[cidi_var]])]
  )

  # As in the original do-file, age contributes to the person-level mean
  # only in waves with a non-missing CIDI probability.
  CIDI_wide[[age_var]][is.na(CIDI_wide[[prob_var]])] <- NA_real_
}

probability_vars <- paste0("Depr_prob_", cidi_probability_years)
age_vars <- paste0("age_", cidi_probability_years)

CIDI_analysis <- CIDI_wide %>%
  mutate(
    Depr_prob = rowMeans(across(all_of(probability_vars)), na.rm = TRUE),
    age = rowMeans(across(all_of(age_vars)), na.rm = TRUE),
    Depr_prob = replace(Depr_prob, is.nan(Depr_prob), NA_real_),
    age = replace(age, is.nan(age), NA_real_),
    age2 = age^2,
    male_age = male * age,
    male_age2 = male * age2
  )

# Convert each wave's CIDI score to depression probability.
cidi_probability <- c(
  `-1` = 0,
  `0` = 0.0001,
  `1` = 0.0568,
  `2` = 0.2352,
  `3` = 0.5542,
  `4` = 0.8125,
  `5` = 0.8895,
  `6` = 0.8895,
  `7` = 0.9083
)

CIDI <- CIDI_long %>%
  mutate(
    pheno = unname(cidi_probability[as.character(CIDI)]),
    age2 = age^2,
    male_age = male * age,
    male_age2 = male * age2
  ) %>%
  select(HHID, PN, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Residualise within wave, standardise within wave, and average across waves.
residualise.average.save(
  data = CIDI,
  average = TRUE,
  age_residualise = TRUE,
  sexres = TRUE,
  rand = FALSE,
  name = "DEP"
)


########################################################
############# Construct phenotype: activity ############
########################################################

# select and rename variables
activity_2004 <- data_2004 %>%
 select(
   HHID, PN,
   pheno1_2004=JC223,
   pheno2_2004=JC224,
   pheno3_2004=JC225,
   sex, mob, yob, moi_2004, yoi_2004
 )

activity_2006 <- data_2006 %>%
 select(
   HHID, PN,
   pheno1_2006=KC223,
   pheno2_2006=KC224,
   pheno3_2006=KC225,
   sex, mob, yob, moi_2006, yoi_2006
 )

activity_2008 <- data_2008 %>%
 select(
   HHID, PN,
   pheno1_2008=LC223,
   pheno2_2008=LC224,
   pheno3_2008=LC225,
   sex, mob, yob, moi_2008, yoi_2008
 )

activity_2010 <- data_2010 %>%
 select(
   HHID, PN,
   pheno1_2010=MC223,
   pheno2_2010=MC224,
   pheno3_2010=MC225,
   sex, mob, yob, moi_2010, yoi_2010
 )

activity_2012 <- data_2012 %>%
 select(
   HHID, PN,
   pheno1_2012=NC223,
   pheno2_2012=NC224,
   pheno3_2012=NC225,
   sex, mob, yob, moi_2012, yoi_2012
 )

activity_2014 <- data_2014 %>%
 select(
   HHID, PN,
   pheno1_2014=OC223,
   pheno2_2014=OC224,
   pheno3_2014=OC225,
   sex, mob, yob, moi_2014, yoi_2014
 )

activity_2016 <- data_2016 %>%
  select(
    HHID, PN,
    pheno1_2016=PC223,
    pheno2_2016=PC224,
    pheno3_2016=PC225,
    sex, mob, yob, moi_2016, yoi_2016
  )

activity_2018 <- data_2018 %>%
  select(
    HHID, PN,
    pheno1_2018=QC223,
    pheno2_2018=QC224,
    pheno3_2018=QC225,
    sex, mob, yob, moi_2018, yoi_2018
  )

# transform data for residualising
activity <- activity_2004 %>%
 full_join(activity_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 select(HHID, PN,
        pheno1_2004, pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016, pheno1_2018,
        pheno2_2004, pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016, pheno2_2018,
        pheno3_2004, pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016, pheno3_2018,
        sex,
        mob,
        yob,
        moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
        yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
 gather(key="wave", value="value",
        paste0("pheno1_", seq(2004, 2018, 2)),
        paste0("pheno2_", seq(2004, 2018, 2)),
        paste0("pheno3_", seq(2004, 2018, 2)),
        paste0("moi_",    seq(2004, 2018, 2)),
        paste0("yoi_",    seq(2004, 2018, 2))) %>%
 separate("wave", c("var", "col")) %>%
 spread("var", "value") %>%
 gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:3)) %>%
 spread("pheno_number", "pheno_value") %>%
 mutate(
   pheno1 = case_when(
     pheno1 == 4 ~ 0, # hardly ever or never
     pheno1 == 3 ~ 0.5, # 1-3 per month, or ~0.5 per week
     pheno1 == 2 ~ 1, # weekly
     pheno1 == 1 ~ 3, # >1 per week
     pheno1 == 7 ~ 7 # daily
   ),
   pheno2 = case_when(
     pheno2 == 4 ~ 0, # hardly ever or never
     pheno2 == 3 ~ 0.5, # 1-3 per month, or ~0.5 per week
     pheno2 == 2 ~ 1, # weekly
     pheno2 == 1 ~ 3, # >1 per week
     pheno2 == 7 ~ 7 # daily
   ),
   pheno3 = case_when(
     pheno3 == 4 ~ 0, # hardly ever or never
     pheno3 == 3 ~ 0.5, # 1-3 per month, or ~0.5 per week
     pheno3 == 2 ~ 1, # weekly
     pheno3 == 1 ~ 3, # >1 per week
     pheno3 == 7 ~ 7 # daily
   ),
   pheno = (8 * pheno1) + (4 * pheno2) + (2 * pheno1)
 ) %>%
 arrange(HHID, PN, col, pheno) %>%
 group_by(HHID, PN, col) %>%
 mutate(dob = yob + (mob/12),
        doi = yoi + (moi/12),
        age = doi - dob,
        age2 = age^2,
        male = 2 - sex,
        male_age = male * age,
        male_age2 = male * age2) %>%
 ungroup() %>%
 select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
 drop_na()

# residualise, average, save
residualise.average.save(data=activity, average=T, age_residualise=T, name="ACTIVITY")


########################################################
############### Construct phenotype: adhd ##############
########################################################

# select and rename relevant variables
adhd_2016 <- data_2016 %>%
  select(HHID, PN,
         sex, mob, yob, moi_2016, yoi_2016,
         contains("PV0")) %>%
  rename(moi = moi_2016,
         yoi = yoi_2016,
         pheno1 = PV001,
         pheno2 = PV002,
         pheno3 = PV003,
         pheno4 = PV004,
         pheno5 = PV005,
         pheno6 = PV006,
         pheno6 = PV006,
         pheno7 = PV007,
         pheno8 = PV008,
         pheno9 = PV009,
         pheno10 = PV010,
         pheno11 = PV011,
         pheno12 = PV012,
         pheno13 = PV013,
         pheno14 = PV014,
         pheno15 = PV015,
         pheno16 = PV016,
         pheno17 = PV017,
         pheno18 = PV018) %>%
  select(-contains("PV"))

# transform data for residualising
adhd <- adhd_2016 %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:18)) %>%
  arrange(HHID, PN, pheno_number) %>%
  group_by(HHID, PN) %>%
  mutate(
    pheno_value = case_when(pheno_value == 1 ~ 1,
                            pheno_value == 5 ~ 0),
    pheno_count = sum(!is.na(pheno_value)),
    pheno = case_when(pheno_count >= 10 ~ mean(pheno_value, na.rm = T)),
    dob = yob + (mob/12),
    doi = yoi + (moi/12),
    age = doi - dob,
    age2 = age^2,
    male = 2 - sex,
    male_age = male * age,
    male_age2 = male * age2
  ) %>%
  ungroup() %>%
  group_by(HHID, PN, age, age2, male, male_age, male_age2) %>%
  summarise(
    pheno = mean(pheno),
    pheno_count = mean(pheno_count)
  ) %>%
  ungroup() %>%
  select(HHID, PN, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
adhd <- adhd %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- adhd %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in adhd:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(adhd)) 

# residualise, average, save
residualise.average.save(data=adhd, average=F, age_residualise=T, name="ADHD")


########################################################
########### Construct phenotype: adventurous ###########
########################################################

# select and rename variables
adventurous_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KLB033Y,
         sex, mob, yob, moi_2006, yoi_2006)

adventurous_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LLB033Y,
         sex, mob, yob, moi_2008, yoi_2008)

adventurous_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MLB033Z_4,
         sex, mob, yob, moi_2010, yoi_2010)

adventurous_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NLB033Z_4,
         sex, mob, yob, moi_2012, yoi_2012)

adventurous_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OLB031Z_4,
         sex, mob, yob, moi_2014, yoi_2014)

adventurous_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PLB031Z_4,
         sex, mob, yob, moi_2016, yoi_2016)

# transform data for residualising
adventurous <- adventurous_2006 %>%
  full_join(adventurous_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(adventurous_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(adventurous_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(adventurous_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(adventurous_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2006, 2016, 2)),
         paste0("moi_",   seq(2006, 2016, 2)),
         paste0("yoi_",   seq(2006, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==1 ~ 3, pheno==2 ~ 2, pheno==3 ~ 1, pheno==4 ~ 0)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=adventurous, average=T, age_residualise=T, name="ADVENTURE")

##################################################################
######## Construct phenotype: age at smoking initiation ##########
#################################################################

attr(data_1998$F1271, "label")
summary(data_1998$F1271)
table(data_1998$F1271)
glimpse(data_1998$F1271)
print(freq_table <- table(data_1998$F1271))

# select and rename relevant variables
ASI_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=F1271,
         sex, mob, yob, moi_1998, yoi_1998)

ASI_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G1404,
         sex, mob, yob, moi_2000, yoi_2000)

ASI_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC120,
         sex, mob, yob, moi_2002, yoi_2002)

ASI_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC120,
         sex, mob, yob, moi_2004, yoi_2004)

ASI_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KC120,
         sex, mob, yob, moi_2006, yoi_2006)

ASI_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LC120,
         sex, mob, yob, moi_2008, yoi_2008)

ASI_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC120,
         sex, mob, yob, moi_2010, yoi_2010)

ASI_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC120,
         sex, mob, yob, moi_2012, yoi_2012)

ASI_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC120,
         sex, mob, yob, moi_2014, yoi_2014)

ASI_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC120,
         sex, mob, yob, moi_2016, yoi_2016)

ASI_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC120,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
ASI <- ASI_1998 %>%
  full_join(ASI_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(ASI_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno_1998, pheno_2000, pheno_2002, pheno_2004, pheno_2006,
         pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,

         sex, mob, yob,

         moi_1998, moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,

         yoi_1998, yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(1998, 2018, 2)),
         paste0("moi_",   seq(1998, 2018, 2)),
         paste0("yoi_",   seq(1998, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = ifelse(!is.na(pheno) & pheno < 98, max(pheno, na.rm = TRUE), NA),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
ASI <- ASI %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- ASI %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in ASI:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(ASI)) 

# residualise, average, save
residualise.average.save(data=ASI, average=F, age_residualise=F, name="ASI")

########################################################
########## Construct phenotype: agreeableness ##########
########################################################

# select and rename variables
agree_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno1_2006=KLB033B, pheno2_2006=KLB033F, pheno3_2006=KLB033J, pheno4_2006=KLB033O, pheno5_2006=KLB033V,
         sex, mob, yob, moi_2006, yoi_2006)

agree_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno1_2008=LLB033B, pheno2_2008=LLB033F, pheno3_2008=LLB033J, pheno4_2008=LLB033O, pheno5_2008=LLB033V,
         sex, mob, yob, moi_2008, yoi_2008)

agree_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno1_2010=MLB033B, pheno2_2010=MLB033G, pheno3_2010=MLB033K, pheno4_2010=MLB033P, pheno5_2010=MLB033Y,
         sex, mob, yob, moi_2010, yoi_2010)

agree_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno1_2012=NLB033B, pheno2_2012=NLB033G, pheno3_2012=NLB033K, pheno4_2012=NLB033P, pheno5_2012=NLB033Y,
         sex, mob, yob, moi_2012, yoi_2012)

agree_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno1_2014=OLB031B, pheno2_2014=OLB031G, pheno3_2014=OLB031K, pheno4_2014=OLB031P, pheno5_2014=OLB031Y,
         sex, mob, yob, moi_2014, yoi_2014)

agree_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno1_2016=PLB031B, pheno2_2016=PLB031G, pheno3_2016=PLB031K, pheno4_2016=PLB031P, pheno5_2016=PLB031Y,
         sex, mob, yob, moi_2016, yoi_2016)

# transform data for residualising
agree <- agree_2006 %>%
  full_join(agree_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(agree_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(agree_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(agree_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(agree_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016,
         pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016,
         pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016,
         pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016,
         pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016,
         sex,
         mob,
         yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno1_", seq(2006, 2016, 2)),
         paste0("pheno2_", seq(2006, 2016, 2)),
         paste0("pheno3_", seq(2006, 2016, 2)),
         paste0("pheno4_", seq(2006, 2016, 2)),
         paste0("pheno5_", seq(2006, 2016, 2)),
         paste0("moi_",    seq(2006, 2016, 2)),
         paste0("yoi_",    seq(2006, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:5)) %>%
  arrange(HHID, PN, col, pheno_number) %>%
  group_by(HHID, PN, col) %>%
  mutate(pheno_count = 5 - sum(is.na(pheno_value)),
         pheno = case_when(pheno_count > 2 ~ mean(pheno_value, na.rm = T)),
         dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>%
  ungroup() %>%
  group_by(HHID, PN, col, age, age2, male, male_age, male_age2) %>%
  summarise(pheno = 5 - mean(pheno), # reverse code
            pheno_count = mean(pheno_count)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
#residualise.average.save(data=agree, average=T, age_residualise=T, name="AGREE")


########################################################
############## Construct phenotype: asthma #############
########################################################

# select and rename relevant variables
asthma_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LB105,
         sex, mob, yob, moi_2008, yoi_2008)

asthma_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MB105,
         sex, mob, yob, moi_2010, yoi_2010)

asthma_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NB105,
         sex, mob, yob, moi_2012, yoi_2012)

asthma_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OB105,
         sex, mob, yob, moi_2014, yoi_2014)

asthma_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PB105,
         sex, mob, yob, moi_2016, yoi_2016)

asthma_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QB105,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
asthma <- asthma_2008 %>%
  full_join(asthma_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,
         sex,
         mob,
         yob,
         moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2008, 2018, 2)),
         paste0("moi_",   seq(2008, 2018, 2)),
         paste0("yoi_",   seq(2008, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
asthma <- asthma %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- asthma %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in asthma:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(asthma)) 

# residualise, average, save
residualise.average.save(data=asthma, average=F, age_residualise=F, name="ASTHMA")


########################################################
############# Construct phenotype: hayfever ############
########################################################

# # select and rename relevant variables
# hayfever_2008 <- data_2008 %>%
#   select(HHID, PN,
#          pheno_2008=LB107,
#          sex, mob, yob, moi_2008, yoi_2008)

# hayfever_2010 <- data_2010 %>%
#   select(HHID, PN,
#          pheno_2010=MB107,
#          sex, mob, yob, moi_2010, yoi_2010)

# hayfever_2012 <- data_2012 %>%
#   select(HHID, PN,
#          pheno_2012=NB107,
#          sex, mob, yob, moi_2012, yoi_2012)

# hayfever_2014 <- data_2014 %>%
#   select(HHID, PN,
#          pheno_2014=OB107,
#          sex, mob, yob, moi_2014, yoi_2014)

# hayfever_2016 <- data_2016 %>%
#   select(HHID, PN,
#          pheno_2016=PB107,
#          sex, mob, yob, moi_2016, yoi_2016)

# hayfever_2018 <- data_2018 %>%
#   select(HHID, PN,
#          pheno_2018=QB107,
#          sex, mob, yob, moi_2018, yoi_2018)

# # transform data for residualising
# hayfever <- hayfever_2008 %>%
#   full_join(hayfever_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#   full_join(hayfever_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#   full_join(hayfever_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#   full_join(hayfever_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#   full_join(hayfever_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#   select(HHID, PN,
#          pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,
#          sex,
#          mob,
#          yob,
#          moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
#          yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
#   gather(key="wave", value="value",
#          paste0("pheno_", seq(2008, 2018, 2)),
#          paste0("moi_",   seq(2008, 2018, 2)),
#          paste0("yoi_",   seq(2008, 2018, 2))) %>%
#   separate("wave", c("var", "col")) %>%
#   spread("var", "value") %>%
#   mutate(dob = yob + (mob/12),
#          dob2 = dob^2,
#          male = 2 - sex,
#          male_dob = male * dob,
#          male_dob2 = male * dob2,
#          pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
#   group_by(HHID, PN) %>%
#   mutate(pheno = max(pheno, na.rm=T),
#          rn = row_number()) %>%
#   ungroup() %>%
#   filter(rn==1 & pheno!="-Inf") %>%
#   select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
#   drop_na()

# # residualise, average, save
# residualise.average.save(data=hayfever, average=F, age_residualise=F, name="HAYFEVER")


########################################################
##### Construct phenotype: asthma-(eczema-)hayfever ####
########################################################

# get asthmaeczrhi from asthma and hayfever above
asthma_hayfever <- asthma %>%
  rename(asthma = pheno) %>%
  full_join(hayfever, by=c("HHID", "PN",
                           "dob", "dob2", "male", "male_dob", "male_dob2")) %>%
  rename(hayfever = pheno) %>%
  drop_na(asthma, hayfever) %>%
  mutate(pheno = case_when(asthma + hayfever > 0 ~ 1,
                           asthma + hayfever == 0 ~ 0)) %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# residualise, average, save
# residualise.average.save(data=asthma_hayfever, average=F, age_residualise=F, name="ASTECZRHI")


########################################################
############## Construct phenotype: audit ##############
########################################################

# HRS Qs:
#
#    "In your entire life, have you had at least 12 drinks of any type of alcoholic beverage?" AUDIT 0 if no
#    "Do you ever drink any alcoholic beverages such as beer, wine, or liquor?" 0 if no
#1    "In the last three months, on average, how many days per week have you had any alcohol to drink?" 0 if [0,1], 1 if [>1].
#2    "In the last three months, on the days you drink, about how many drinks do you have?" 0 if [0-4], 1 if [>4].
#3    "In the last three months, on how many days have you had four or more drinks on one occasion?" 0 if [0-3], 1 if [>3].
#    Rest are 0-1. Use


# "Have you ever taken a drink first thing in the morning to steady your nerves or get rid of a hangover?"
# "Have you ever felt bad or guilty about drinking?"
# "Have you ever felt that you should cut down on drinking?"
# "Have people ever annoyed you by criticizing your drinking?".

# select and rename variables
audit_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno1_2002=HC134,
         pheno2_2002=HC128,

         pheno3_2002=HC129,
         pheno4_2002=HC130,
         pheno5_2002=HC131,
         pheno6_2002=HC138,
         pheno7_2002=HC137,
         pheno8_2002=HC135,
         pheno9_2002=HC136,
         sex, mob, yob, moi_2002, yoi_2002)

audit_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno1_2004=JC134,
         pheno2_2004=JC128,

         pheno3_2004=JC129,
         pheno4_2004=JC130,
         pheno5_2004=JC131,
         pheno6_2004=JC138,
         pheno7_2004=JC137,
         pheno8_2004=JC135,
         pheno9_2004=JC136,
         sex, mob, yob, moi_2004, yoi_2004)

audit_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno1_2006=KC134,
         pheno2_2006=KC128,

         pheno3_2006=KC129,
         pheno4_2006=KC130,
         pheno5_2006=KC131,
         pheno6_2006=KC138,
         pheno7_2006=KC137,
         pheno8_2006=KC135,
         pheno9_2006=KC136,
         sex, mob, yob, moi_2006, yoi_2006)

audit_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno1_2008=LC134,
         pheno2_2008=LC128,

         pheno3_2008=LC129,
         pheno4_2008=LC130,
         pheno5_2008=LC131,
         pheno6_2008=LC138,
         pheno7_2008=LC137,
         pheno8_2008=LC135,
         pheno9_2008=LC136,
         sex, mob, yob, moi_2008, yoi_2008)

audit_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno1_2010=MC134,
         pheno2_2010=MC128,

         pheno3_2010=MC129,
         pheno4_2010=MC130,
         pheno5_2010=MC131,
         pheno6_2010=MC138,
         pheno7_2010=MC137,
         pheno8_2010=MC135,
         pheno9_2010=MC136,
         sex, mob, yob, moi_2010, yoi_2010)

audit_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno1_2012=NC134,
         pheno2_2012=NC128,

         pheno3_2012=NC129,
         pheno4_2012=NC130,
         pheno5_2012=NC131,
         pheno6_2012=NC138,
         pheno7_2012=NC137,
         pheno8_2012=NC135,
         pheno9_2012=NC136,
         sex, mob, yob, moi_2012, yoi_2012)

audit_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno1_2014=OC134,
         pheno2_2014=OC128,

         pheno3_2014=OC129,
         pheno4_2014=OC130,
         pheno5_2014=OC131,
         pheno6_2014=OC138,
         pheno7_2014=OC137,
         pheno8_2014=OC135,
         pheno9_2014=OC136,
         sex, mob, yob, moi_2014, yoi_2014)

audit_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno1_2016=PC134,
         pheno2_2016=PC128,

         pheno3_2016=PC129,
         pheno4_2016=PC130,
         pheno5_2016=PC131,
         pheno6_2016=PC138,
         pheno7_2016=PC137,
         pheno8_2016=PC135,
         pheno9_2016=PC136,
         sex, mob, yob, moi_2016, yoi_2016)

audit_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno1_2018=QC134,
         pheno2_2018=QC128,

         pheno3_2018=QC129,
         pheno4_2018=QC130,
         pheno5_2018=QC131,
         pheno6_2018=QC138,
         pheno7_2018=QC137,
         pheno8_2018=QC135,
         pheno9_2018=QC136,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
audit <- audit_2002 %>%
  full_join(audit_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno1_2002, pheno1_2004, pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016, pheno1_2018,
         pheno2_2002, pheno2_2004, pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016, pheno2_2018,
         pheno3_2002, pheno3_2004, pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016, pheno3_2018,
         pheno4_2002, pheno4_2004, pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016, pheno4_2018,
         pheno5_2002, pheno5_2004, pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016, pheno5_2018,
         pheno6_2002, pheno6_2004, pheno6_2006, pheno6_2008, pheno6_2010, pheno6_2012, pheno6_2014, pheno6_2016, pheno6_2018,
         pheno7_2002, pheno7_2004, pheno7_2006, pheno7_2008, pheno7_2010, pheno7_2012, pheno7_2014, pheno7_2016, pheno7_2018,
         pheno8_2002, pheno8_2004, pheno8_2006, pheno8_2008, pheno8_2010, pheno8_2012, pheno8_2014, pheno8_2016, pheno8_2018,
         pheno9_2002, pheno9_2004, pheno9_2006, pheno9_2008, pheno9_2010, pheno9_2012, pheno9_2014, pheno9_2016, pheno9_2018,
         sex,
         mob,
         yob,
         moi_2002, moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2002, yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno1_", seq(2002, 2018, 2)),
         paste0("pheno2_", seq(2002, 2018, 2)),
         paste0("pheno3_", seq(2002, 2018, 2)),
         paste0("pheno4_", seq(2002, 2018, 2)),
         paste0("pheno5_", seq(2002, 2018, 2)),
         paste0("pheno6_", seq(2002, 2018, 2)),
         paste0("pheno7_", seq(2002, 2018, 2)),
         paste0("pheno8_", seq(2002, 2018, 2)),
         paste0("pheno9_", seq(2002, 2018, 2)),
         paste0("moi_",    seq(2002, 2018, 2)),
         paste0("yoi_",    seq(2002, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:9)) %>%
  spread("pheno_number", "pheno_value") %>%
  mutate(
    pheno1 = case_when(
      pheno1 %in% 1:2 ~ 1,
      pheno1 %in% 5:6 ~ 0
    ),
    pheno2 = case_when(
      pheno2 %in% 1:2 ~ 1,
      pheno2 %in% 3:6 ~ 0
    ),
    # "In the last three months, on average, how many days per week have you had any alcohol to drink?"
    # 0 if [0,1], 1 if [>1].
    pheno3 = case_when(
      #pheno3 %in% 0:7 ~ as.numeric(pheno3),
      pheno3 %in% 0:1 ~ 0,
      pheno3 %in% 2:7 ~ 1
    ),
    # "In the last three months, on the days you drink, about how many drinks do you have?"
    # 0 if [0/1-4], 1 if [>4].
    pheno4 = case_when(
      #pheno4 %in% 0:50 ~ as.numeric(pheno4),
      pheno4 %in% 0:4 ~ 0,
      pheno4 %in% 5:50 ~ 1
    ),
    # "In the last three months, on how many days have you had four or more drinks on one occasion?"
    # 0 if [0-3], 1 if [>3].
    pheno5 = case_when(
      #pheno5 %in% 0:97 ~ as.numeric(pheno5),
      pheno5 %in% 0:3 ~ 0,
      pheno5 %in% 4:97 ~ 1
    ),
    pheno6 = case_when(
      pheno6 %in% 1:2 ~ 1,
      pheno6 %in% 5:6 ~ 0
    ),
    pheno7 = case_when(
      pheno7 %in% 1:2 ~ 1,
      pheno7 %in% 5:6 ~ 0
    ),
    pheno8 = case_when(
      pheno8 %in% 1:2 ~ 1,
      pheno8 %in% 5:6 ~ 0
    ),
    pheno9 = case_when(
      pheno9 %in% 1:2 ~ 1,
      pheno9 %in% 5:6 ~ 0
    )
  ) %>%
  rowwise() %>%
  mutate(
    pheno_NA = sum(is.na(c(pheno3, pheno4, pheno5, pheno6, pheno7, pheno8, pheno9))),
    pheno_filter = case_when(
      # explicitly code all cases *except when both values are NA* to avoid false-0s
      pheno1 == 0 ~ 0,
      pheno2 == 0 ~ 0,
      pheno1 > 0 ~ 1,
      pheno2 > 0 ~ 1
    ),
    pheno_sum = sum(c(pheno3, pheno4, pheno5, pheno6, pheno7, pheno8, pheno9), na.rm = F)
  ) %>%
  ungroup() %>%
  mutate(
  #  pheno_count = 7 - pheno_NA,
    pheno = case_when(
      pheno_filter == 1 ~ pheno_sum,
      pheno_filter == 0 ~ 0
    )
  ) %>%
  #tail(select(audit, HHID, PN, pheno_filter, pheno3, pheno4, pheno5, pheno_sum, pheno_count, pheno), n=20)
  arrange(HHID, PN, col, pheno) %>%
  group_by(HHID, PN, col) %>%
  mutate(
    dob = yob + (mob/12),
    doi = yoi + (moi/12),
    age = doi - dob,
    age2 = age^2,
    male = 2 - sex,
    male_age = male * age,
    male_age2 = male * age2
  ) %>%
  ungroup() %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=audit, average=T, age_residualise=T, name="AUDIT")

##########################################################
############## Construct phenotype: diabetes #############
##########################################################

# select and rename relevant variables
diabetes_1995 <- data_1995 %>%
  select(HHID, PN,
         pheno_1995=D788,
         sex, mob, yob, moi_1995, yoi_1995)

diabetes_1996 <- data_1996 %>%
  select(HHID, PN,
         pheno_1996=E788,
         sex, mob, yob, moi_1996, yoi_1996)

diabetes_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=F1116,
         sex, mob, yob, moi_1998, yoi_1998)

diabetes_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G1245,
         sex, mob, yob, moi_2000, yoi_2000)

diabetes_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC010,
         sex, mob, yob, moi_2002, yoi_2002)

diabetes_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC010,
         sex, mob, yob, moi_2004, yoi_2004)

diabetes_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KC010,
         sex, mob, yob, moi_2006, yoi_2006)

diabetes_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LC010,
         sex, mob, yob, moi_2008, yoi_2008)

diabetes_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC010,
         sex, mob, yob, moi_2010, yoi_2010)

diabetes_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC010,
         sex, mob, yob, moi_2012, yoi_2012)

diabetes_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC010,
         sex, mob, yob, moi_2014, yoi_2014)

diabetes_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC010,
         sex, mob, yob, moi_2016, yoi_2016)

diabetes_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC010,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
diabetes <- diabetes_1995 %>%
  full_join(diabetes_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(diabetes_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,

         sex,
         mob,
         yob,

         moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,

         yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1995, seq(1996, 2018, 2))),
         paste0("moi_",   c(1995, seq(1996, 2018, 2))),
         paste0("yoi_",   c(1995, seq(1996, 2018, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
diabetes <- diabetes %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- diabetes %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in diabetes:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(diabetes)) 

# residualise, average, save
residualise.average.save(data=diabetes, average=F, age_residualise=F, name="T2D")


#############################################################
######### Construct phenotype: bipolar disorder #############
#############################################################

# select and rename relevant variables
bipolar_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OV353,
         sex, mob, yob, moi_2014, yoi_2014)

# transform data for residualising
bipolar <- bipolar_2014 %>%
  select(HHID, PN, pheno_2014, sex, mob, yob, moi_2014, yoi_2014) %>%
  gather(key="wave", value="value",
         paste0("pheno_", 2014),
         paste0("moi_",   2014),
         paste0("yoi_",   2014)) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = ifelse(!is.na(pheno) & pheno < 8, max(pheno, na.rm = TRUE), NA),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
bipolar <- bipolar %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- bipolar %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in bipolar:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(bipolar)) 

# residualise, average, save
residualise.average.save(data=bipolar, average=F, age_residualise=F, name="BIPOLAR")

#####################################################################
######## Construct phenotype: Blood pressure - Diastolic ############
#####################################################################

# 2006 Wave
Diastolic_2006 <- data_2006 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(KI860 > 0 & KI860 < 900, KI860, NA),
         measure2 = ifelse(KI865 > 0 & KI865 < 900, KI865, NA),
         measure3 = ifelse(KI870 > 0 & KI870 < 900, KI870, NA),
         pheno_2006 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2006, sex, mob, yob, moi_2006, yoi_2006)

# 2008 Wave
Diastolic_2008 <- data_2008 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(LI860 > 0 & LI860 < 900, LI860, NA),
         measure2 = ifelse(LI865 > 0 & LI865 < 900, LI865, NA),
         measure3 = ifelse(LI870 > 0 & LI870 < 900, LI870, NA),
         pheno_2008 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2008, sex, mob, yob, moi_2008, yoi_2008)

# 2010 Wave
Diastolic_2010 <- data_2010 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(MI860 > 0 & MI860 < 900, MI860, NA),
         measure2 = ifelse(MI865 > 0 & MI865 < 900, MI865, NA),
         measure3 = ifelse(MI870 > 0 & MI870 < 900, MI870, NA),
         pheno_2010 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2010, sex, mob, yob, moi_2010, yoi_2010)

# 2012 Wave
Diastolic_2012 <- data_2012 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(NI860 > 0 & NI860 < 900, NI860, NA),
         measure2 = ifelse(NI865 > 0 & NI865 < 900, NI865, NA),
         measure3 = ifelse(NI870 > 0 & NI870 < 900, NI870, NA),
         pheno_2012 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2012, sex, mob, yob, moi_2012, yoi_2012)

# 2014 Wave
Diastolic_2014 <- data_2014 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(OI860 > 0 & OI860 < 900, OI860, NA),
         measure2 = ifelse(OI865 > 0 & OI865 < 900, OI865, NA),
         measure3 = ifelse(OI870 > 0 & OI870 < 900, OI870, NA),
         pheno_2014 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2014, sex, mob, yob, moi_2014, yoi_2014)

# 2018 Wave
Diastolic_2018 <- data_2018 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(QI860 > 0 & QI860 < 900, QI860, NA),
         measure2 = ifelse(QI865 > 0 & QI865 < 900, QI865, NA),
         measure3 = ifelse(QI870 > 0 & QI870 < 900, QI870, NA),
         pheno_2018 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2018, sex, mob, yob, moi_2018, yoi_2018)

# --- Combine and Transform Diastolic Data for Residualising ---
Diastolic <- Diastolic_2006 %>%
  full_join(Diastolic_2008, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Diastolic_2010, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Diastolic_2012, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Diastolic_2014, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Diastolic_2018, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2018,
         sex, mob, yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2018,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2018) %>%
  gather(key = "wave", value = "value",
         paste0("pheno_", c(seq(2006, 2014, 2), 2018)),
         paste0("moi_",   c(seq(2006, 2014, 2), 2018)),
         paste0("yoi_",   c(seq(2006, 2014, 2), 2018))) %>%
  separate("wave", c("var", "col"), sep = "_") %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob / 12),
         doi = yoi + (moi / 12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>%
  select(HHID, PN, wave = col, pheno, age, age2, male, male_age, male_age2, yoi, moi) %>%
  drop_na() %>%
  # Exclude values outside the plausible range for diastolic BP (in mm Hg):
  filter(pheno >= 50 & pheno <= 150)

# residualise, average, save
residualise.average.save(data=Diastolic, average=T, age_residualise=T, name="BPdia")

#####################################################################
######## Construct phenotype: Blood pressure - Systolic #############
#####################################################################

# 2006 Wave
Systolic_2006 <- data_2006 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(KI859 > 0 & KI859 < 900, KI859, NA),
         measure2 = ifelse(KI864 > 0 & KI864 < 900, KI864, NA),
         measure3 = ifelse(KI869 > 0 & KI869 < 900, KI869, NA),
         pheno_2006 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2006, sex, mob, yob, moi_2006, yoi_2006)

# 2008 Wave
Systolic_2008 <- data_2008 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(LI859 > 0 & LI859 < 900, LI859, NA),
         measure2 = ifelse(LI864 > 0 & LI864 < 900, LI864, NA),
         measure3 = ifelse(LI869 > 0 & LI869 < 900, LI869, NA),
         pheno_2008 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2008, sex, mob, yob, moi_2008, yoi_2008)

# 2010 Wave
Systolic_2010 <- data_2010 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(MI859 > 0 & MI859 < 900, MI859, NA),
         measure2 = ifelse(MI864 > 0 & MI864 < 900, MI864, NA),
         measure3 = ifelse(MI869 > 0 & MI869 < 900, MI869, NA),
         pheno_2010 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2010, sex, mob, yob, moi_2010, yoi_2010)

# 2012 Wave
Systolic_2012 <- data_2012 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(NI859 > 0 & NI859 < 900, NI859, NA),
         measure2 = ifelse(NI864 > 0 & NI864 < 900, NI864, NA),
         measure3 = ifelse(NI869 > 0 & NI869 < 900, NI869, NA),
         pheno_2012 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2012, sex, mob, yob, moi_2012, yoi_2012)

# 2014 Wave
Systolic_2014 <- data_2014 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(OI859 > 0 & OI859 < 900, OI859, NA),
         measure2 = ifelse(OI864 > 0 & OI864 < 900, OI864, NA),
         measure3 = ifelse(OI869 > 0 & OI869 < 900, OI869, NA),
         pheno_2014 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2014, sex, mob, yob, moi_2014, yoi_2014)

# 2018 Wave
Systolic_2018 <- data_2018 %>%
  rowwise() %>%
  mutate(measure1 = ifelse(QI859 > 0 & QI859 < 900, QI859, NA),
         measure2 = ifelse(QI864 > 0 & QI864 < 900, QI864, NA),
         measure3 = ifelse(QI869 > 0 & QI869 < 900, QI869, NA),
         pheno_2018 = mean(c(measure1, measure2, measure3), na.rm = TRUE)) %>%
  select(HHID, PN, pheno_2018, sex, mob, yob, moi_2018, yoi_2018)

# --- Transform Data for Residualising ---
Systolic <- Systolic_2006 %>%
  full_join(Systolic_2008, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Systolic_2010, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Systolic_2012, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Systolic_2014, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Systolic_2018, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2018,
         sex, mob, yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2018,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2018) %>%
  gather(key = "wave", value = "value",
         paste0("pheno_", c(seq(2006, 2014, 2), 2018)),
         paste0("moi_",   c(seq(2006, 2014, 2), 2018)),
         paste0("yoi_",   c(seq(2006, 2014, 2), 2018))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob / 12),
         doi = yoi + (moi / 12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>%
  select(HHID, PN, wave = col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na() %>%
  # Exclude values outside the plausible range for systolic BP (in mm Hg):
  filter(pheno >= 70 & pheno <= 270)
  
# Residualise, average, save
residualise.average.save(data=Systolic, average=T, age_residualise=T, name="BPsys")

#####################################################################
######## Construct phenotype: Blood pressure - Pulse #############
#####################################################################

# Select, average, and rename relevant variables
Pulse_2006 <- data_2006 %>%
  mutate(pheno_2006 = average_three(KI861, KI866, KI871)) %>%
  select(HHID, PN, pheno_2006, sex, mob, yob, moi_2006, yoi_2006)

Pulse_2008 <- data_2008 %>%
  mutate(pheno_2008 = average_three(LI861, LI866, LI871)) %>%
  select(HHID, PN, pheno_2008, sex, mob, yob, moi_2008, yoi_2008)

Pulse_2010 <- data_2010 %>%
  mutate(pheno_2010 = average_three(MI861, MI866, MI871)) %>%
  select(HHID, PN, pheno_2010, sex, mob, yob, moi_2010, yoi_2010)

Pulse_2012 <- data_2012 %>%
  mutate(pheno_2012 = average_three(NI861, NI866, NI871)) %>%
  select(HHID, PN, pheno_2012, sex, mob, yob, moi_2012, yoi_2012)

Pulse_2014 <- data_2014 %>%
  mutate(pheno_2014 = average_three(OI861, OI866, OI871)) %>%
  select(HHID, PN, pheno_2014, sex, mob, yob, moi_2014, yoi_2014)

Pulse_2018 <- data_2018 %>%
  mutate(pheno_2018 = average_three(QI861, QI866, QI871)) %>%
  select(HHID, PN, pheno_2018, sex, mob, yob, moi_2018, yoi_2018)

# Transform data for residualising
Pulse <- Pulse_2006 %>%
  full_join(Pulse_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Pulse_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Pulse_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Pulse_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(Pulse_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2018,
         sex, mob, yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2018,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(seq(2006, 2014, 2), 2018)),
         paste0("moi_",   c(seq(2006, 2014, 2), 2018)),
         paste0("yoi_",   c(seq(2006, 2014, 2), 2018))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob / 12),
         doi = yoi + (moi / 12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na() %>%
  filter(pheno <= 701) # Truncate implausible extreme values (out of range)

# Residualise, average, save
residualise.average.save(data=Pulse, average=T, age_residualise=T, name="PulseRate")

#####################################################################
######## Construct phenotype: Blood pressure - Pulse #############
#####################################################################

# Join Systolic and Diastolic data and compute BPpulse, which is BPsys - BPdia
BPpulse <- Systolic %>%
  inner_join(Diastolic, by=c("HHID", "PN", "wave", "age", "age2", "male", "male_age", "male_age2")) %>%
  mutate(pheno = pheno.x - pheno.y) %>%
  select(HHID, PN, wave, pheno, age, age2, male, male_age, male_age2, yoi, moi) %>%
  drop_na()

# Apply the residualising functions to BPpulse
residualise.average.save(data=BPpulse, average=T, age_residualise=T, name="BPpulse")

####################################################################
########### Construct phenotype: Cancer - breast cancer ############
####################################################################

# Check unique values of 'sex' and their counts to understand the coding
table(data_1992$sex)  # Assuming 'data_1992' can be representative for all years

# Select and rename relevant variables, restrict to females if sex == 2
breastCancer_1992 <- data_1992 %>%
  filter(sex == 2) %>%  # Filtering to include only females
  select(
    HHID, PN,
    pheno_1992 = V340,
    secondCancer_1992 = V347,
    sex, mob, yob, moi_1992, yoi_1992
  )
# Diagnostic 
print(unique(breastCancer_1992$sex))

breastCancer_1993 <- data_1993 %>%
  filter(sex == 2) %>%  # Filtering to include only females
  select(
    HHID, PN,
    pheno_1993 = V231,    
    sex, mob, yob, moi_1993, yoi_1993
  )

breastCancer_1994 <- data_1994 %>%
  filter(sex == 2) %>%  # Filtering to include only females
  select(
    HHID, PN,
    pheno_1994 = W343,   
    sex, mob, yob, moi_1994, yoi_1994
  )

# Transform data for residualising and merge yearly data
breastCancer <- breastCancer_1992 %>%
  full_join(breastCancer_1993, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(breastCancer_1994, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  mutate(
    dob = yob + (mob / 12),
    dob2 = dob^2,
    male = 2 - sex,
    male_dob = male * dob,
    male_dob2 = male * dob2,
    pheno = case_when(
      pheno_1992 == 10 ~ 1, 
      secondCancer_1992 == 10 ~ 1, 
      pheno_1993 == 10 ~ 1, 
      pheno_1994 == 10 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  group_by(HHID, PN) %>%
  mutate(
    pheno = max(pheno, na.rm = TRUE),
    rn = row_number()
  ) %>%
  ungroup() %>%
  filter(rn == 1 & pheno != "-Inf") %>%
  mutate(
    wave = case_when(
      !is.na(pheno_1992) ~ "1992",
      !is.na(pheno_1993) ~ "1993",
      !is.na(pheno_1994) ~ "1994"
    )
  ) %>%
  select(HHID, PN, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
breastCancer <- breastCancer %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- breastCancer %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in breastCancer:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(breastCancer))

# residualise, average, save
residualise.average.save(data=breastCancer, average=F, age_residualise=F, name="BRCA")

####################################################################
########### Construct phenotype: Cancer - prostate cancer ############
####################################################################

# Select and rename relevant variables, restricted to males if sex == 1
prostateCancer_1992 <- data_1992 %>%
  filter(sex == 1) %>%  # Filtering to include only males
  select(
    HHID, PN,
    pheno_1992 = V340,
    secondCancer_1992 = V347,
    sex, mob, yob, moi_1992, yoi_1992
  )
# Diagnostic 
print(unique(prostateCancer_1992$sex))

prostateCancer_1993 <- data_1993 %>%
  filter(sex == 1) %>%  # Filtering to include only males
  select(
    HHID, PN,
    pheno_1993 = V231,    
    sex, mob, yob, moi_1993, yoi_1993
  )

prostateCancer_1994 <- data_1994 %>%
  filter(sex == 1) %>%  # Filtering to include only males
  select(
    HHID, PN,
    pheno_1994 = W343,   
    sex, mob, yob, moi_1994, yoi_1994
  )

# Transform data for residualising and merge yearly data
prostateCancer <- prostateCancer_1992 %>%
  full_join(prostateCancer_1993, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(prostateCancer_1994, by = c("HHID", "PN", "sex", "mob", "yob")) %>%
  mutate(
    dob = yob + (mob / 12),
    dob2 = dob^2,
    male = 2 - sex,  # Calculate the 'male' variable
    male_dob = male * dob,
    male_dob2 = male * dob2,
    pheno = case_when(
      pheno_1992 == 41 ~ 1, 
      secondCancer_1992 == 41 ~ 1, 
      pheno_1993 == 41 ~ 1, 
      pheno_1994 == 41 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  group_by(HHID, PN) %>%
  mutate(
    pheno = max(pheno, na.rm = TRUE),
    rn = row_number()
  ) %>%
  ungroup() %>%
  filter(rn == 1 & pheno != "-Inf") %>%
  mutate(
    wave = case_when(
      !is.na(pheno_1992) ~ "1992",
      !is.na(pheno_1993) ~ "1993",
      !is.na(pheno_1994) ~ "1994"
    )
  ) %>%
  select(HHID, PN, wave, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
prostateCancer <- prostateCancer %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- prostateCancer %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in prostateCancer:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(prostateCancer)) 

# residualise, average, save
residualise.average.save(data=prostateCancer, average=F, age_residualise=F, name="PRCA")
  
########################################################
######## Construct phenotype: conscientiousness ########
########################################################

# select and rename relevant variables
consc_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno1_2006=KLB033D,
         pheno2_2006=KLB033H,
         pheno3_2006=KLB033M,
         pheno4_2006=KLB033T,
         pheno5_2006=KLB033Z,
         sex, mob, yob, moi_2006, yoi_2006)

consc_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno1_2008=LLB033D,
         pheno2_2008=LLB033H,
         pheno3_2008=LLB033M,
         pheno4_2008=LLB033T,
         pheno5_2008=LLB033Z,
         sex, mob, yob, moi_2008, yoi_2008)

consc_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno1_2010=MLB033E,
         pheno2_2010=MLB033I,
         pheno3_2010=MLB033N,
         pheno4_2010=MLB033V,
         pheno5_2010=MLB033Z_5,
         pheno6_2010=MLB033C,
         pheno7_2010=MLB033R,
         pheno8_2010=MLB033X,
         pheno9_2010=MLB033Z,
         pheno10_2010=MLB033Z_6,
         sex, mob, yob, moi_2010, yoi_2010)

consc_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno1_2012=NLB033E,
         pheno2_2012=NLB033I,
         pheno3_2012=NLB033N,
         pheno4_2012=NLB033V,
         pheno5_2012=NLB033Z_5,
         pheno6_2012=NLB033C,
         pheno7_2012=NLB033R,
         pheno8_2012=NLB033X,
         pheno9_2012=NLB033Z,
         pheno10_2012=NLB033Z_6,
         sex, mob, yob, moi_2012, yoi_2012)

consc_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno1_2014=OLB031E,
         pheno2_2014=OLB031I,
         pheno3_2014=OLB031N,
         pheno4_2014=OLB031V,
         pheno5_2014=OLB031Z_5,
         pheno6_2014=OLB031C,
         pheno7_2014=OLB031R,
         pheno8_2014=OLB031X,
         pheno9_2014=OLB031Z_1,
         pheno10_2014=OLB031Z_6,
         sex, mob, yob, moi_2014, yoi_2014)

consc_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno1_2016=PLB031E,
         pheno2_2016=PLB031I,
         pheno3_2016=PLB031N,
         pheno4_2016=PLB031V,
         pheno5_2016=PLB031Z_5,
         pheno6_2016=PLB031C,
         pheno7_2016=PLB031R,
         pheno8_2016=PLB031X,
         pheno9_2016=PLB031Z_1,
         pheno10_2016=PLB031Z_6,
         sex, mob, yob, moi_2016, yoi_2016)

# transform data for residualising
consc <- consc_2006 %>%
  full_join(consc_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(consc_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(consc_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(consc_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(consc_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
        pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016,
        pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016,
        pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016,
        pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016,
        pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016,
        pheno6_2010, pheno6_2012, pheno6_2014, pheno6_2016,
        pheno7_2010, pheno7_2012, pheno7_2014, pheno7_2016,
        pheno8_2010, pheno8_2012, pheno8_2014, pheno8_2016,
        pheno9_2010, pheno9_2012, pheno9_2014, pheno9_2016,
        pheno10_2010, pheno10_2012, pheno10_2014, pheno10_2016,
        sex,
        mob,
        yob,
        moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
        yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
        paste0("pheno1_",  seq(2006, 2016, 2)),
        paste0("pheno2_",  seq(2006, 2016, 2)),
        paste0("pheno3_",  seq(2006, 2016, 2)),
        paste0("pheno4_",  seq(2006, 2016, 2)),
        paste0("pheno5_",  seq(2006, 2016, 2)),
        paste0("pheno6_",  seq(2010, 2016, 2)),
        paste0("pheno7_",  seq(2010, 2016, 2)),
        paste0("pheno8_",  seq(2010, 2016, 2)),
        paste0("pheno9_",  seq(2010, 2016, 2)),
        paste0("pheno10_", seq(2010, 2016, 2)),
        paste0("moi_",     seq(2006, 2016, 2)),
        paste0("yoi_",     seq(2006, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(
    pheno4 = 5 - pheno4,
    pheno6 = 5 - pheno6,
    pheno8 = 5 - pheno8
  ) %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:10)) %>%
  arrange(HHID, PN, col, pheno_number) %>%
  group_by(HHID, PN, col) %>%
  mutate(
    pheno_count = 10 - sum(is.na(pheno_value)),
    pheno = case_when(
      col >  2008 & pheno_count > 4 ~ mean(pheno_value, na.rm = T),
      col <= 2008 & pheno_count > 2 ~ mean(pheno_value, na.rm = T)
    ),
    dob = yob + (mob/12),
    doi = yoi + (moi/12),
    age = doi - dob,
    age2 = age^2,
    male = 2 - sex,
    male_age = male * age,
    male_age2 = male * age2
    ) %>%
  ungroup() %>%
  group_by(HHID, PN, col, age, age2, male, male_age, male_age2) %>%
  summarise(
    pheno = 5 - mean(pheno), # reverse code
    pheno_count = mean(pheno_count)
  )  %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=consc, average=T, age_residualise=T, name="CONSC")


########################################################
############### Construct phenotype: COPD ##############
########################################################

# select and rename variables
copd_1992 <- data_1992 %>%
  select(HHID, PN,
         pheno_1992=V401,
         sex, mob, yob, moi_1992, yoi_1992)

copd_1993 <- data_1993 %>%
  select(HHID, PN,
         pheno_1993=V235,
         sex, mob, yob, moi_1993, yoi_1993)

copd_1994 <- data_1994 %>%
  select(HHID, PN,
         pheno_1994=W362,
         sex, mob, yob, moi_1994, yoi_1994)

copd_1995 <- data_1995 %>%
  select(HHID, PN,
         pheno_1995=D818,
         sex, mob, yob, moi_1995, yoi_1995)

copd_1996 <- data_1996 %>%
  select(HHID, PN,
         pheno_1996=E818,
         sex, mob, yob, moi_1996, yoi_1996)

copd_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=F1146,
         sex, mob, yob, moi_1998, yoi_1998)

copd_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G1279,
         sex, mob, yob, moi_2000, yoi_2000)

copd_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC030,
         sex, mob, yob, moi_2002, yoi_2002)

copd_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC030,
         sex, mob, yob, moi_2004, yoi_2004)

copd_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KC030,
         sex, mob, yob, moi_2006, yoi_2006)

copd_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LC030,
         sex, mob, yob, moi_2008, yoi_2008)

copd_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC030,
         sex, mob, yob, moi_2010, yoi_2010)

copd_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC030,
         sex, mob, yob, moi_2012, yoi_2012)

copd_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC030,
         sex, mob, yob, moi_2014, yoi_2014)

copd_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC030,
         sex, mob, yob, moi_2016, yoi_2016)

copd_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC030,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
copd <- copd_1992 %>%
  full_join(copd_1993, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_1994, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_1995, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(copd_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_1992, pheno_1993, pheno_1994, pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,
         sex,
         mob,
         yob,
         moi_1992, moi_1993, moi_1994, moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016,  moi_2018,
         yoi_1992, yoi_1993, yoi_1994, yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1992:1996, seq(1998, 2018, 2))),
         paste0("moi_",   c(1992:1996, seq(1998, 2018, 2))),
         paste0("yoi_",   c(1992:1996, seq(1998, 2018, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno %in% 1:3 ~ 1, pheno %in% 4:5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
copd <- copd %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- copd %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in copd:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(copd)) 

# residualise, average, save
residualise.average.save(data=copd, average=F, age_residualise=F, name="COPD")


########################################################
############### Construct phenotype: CPD ###############
########################################################

# select and rename relevant variables
CPD_1992 <- data_1992 %>%
  select(HHID, PN,
         pheno_1992=V503,
         sex, mob, yob, moi_1992, yoi_1992)

CPD_1993 <- data_1993 %>%
  select(HHID, PN,
         pheno_1993=V299,
         sex, mob, yob, moi_1993, yoi_1993)

CPD_1994 <- data_1994 %>%
  select(HHID, PN,
         pheno_1994=W453,
         sex, mob, yob, moi_1994, yoi_1994)

CPD_1995 <- data_1995 %>%
  select(HHID, PN,
         pheno_1995=D943,
         sex, mob, yob, moi_1995, yoi_1995)

CPD_1996 <- data_1996 %>%
  select(HHID, PN,
         pheno_1996=E943,
         sex, mob, yob, moi_1996, yoi_1996)

CPD_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=F1268,
         sex, mob, yob, moi_1998, yoi_1998)

CPD_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G1401,
         sex, mob, yob, moi_2000, yoi_2000)

CPD_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC118,
         sex, mob, yob, moi_2002, yoi_2002)

CPD_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC118,
         sex, mob, yob, moi_2004, yoi_2004)

CPD_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KC118,
         sex, mob, yob, moi_2006, yoi_2006)

CPD_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LC118,
         sex, mob, yob, moi_2008, yoi_2008)

CPD_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC118,
         sex, mob, yob, moi_2010, yoi_2010)

CPD_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC118,
         sex, mob, yob, moi_2012, yoi_2012)

CPD_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC118,
         sex, mob, yob, moi_2014, yoi_2014)

CPD_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC118,
         sex, mob, yob, moi_2016, yoi_2016)

CPD_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC118,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
CPD <- CPD_1992 %>%
  full_join(CPD_1993, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_1994, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_1995, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(CPD_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno_1992, pheno_1993, pheno_1994, pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,

         sex, mob, yob,

         moi_1992, moi_1993, moi_1994, moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,

         yoi_1992, yoi_1993, yoi_1994, yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1992:1995, seq(1996, 2018, 2))),
         paste0("moi_",   c(1992:1995, seq(1996, 2018, 2))),
         paste0("yoi_",   c(1992:1995, seq(1996, 2018, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na() %>%
  filter(pheno <= 100) # Truncate implausible extreme values

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
CPD <- CPD %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- CPD %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in CPD:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(CPD)) 

# residualise, average, save
residualise.average.save(data=CPD, average=T, age_residualise=T, name="CPD")


########################################################
########### Construct phenotype: extraversion ##########
########################################################

# select and rename relevant variables
extraversion_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno1_2006=KLB033A,
         pheno2_2006=KLB033E,
         pheno3_2006=KLB033I,
         pheno4_2006=KLB033S,
         pheno5_2006=KLB033W,
         sex, mob, yob, moi_2006, yoi_2006)

extraversion_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno1_2008=LLB033A,
         pheno2_2008=LLB033E,
         pheno3_2008=LLB033I,
         pheno4_2008=LLB033S,
         pheno5_2008=LLB033W,
         sex, mob, yob, moi_2008, yoi_2008)

extraversion_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno1_2010=MLB033A,
         pheno2_2010=MLB033F,
         pheno3_2010=MLB033J,
         pheno4_2010=MLB033U,
         pheno5_2010=MLB033Z_2,
         sex, mob, yob, moi_2010, yoi_2010)

extraversion_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno1_2012=NLB033A,
         pheno2_2012=NLB033F,
         pheno3_2012=NLB033J,
         pheno4_2012=NLB033U,
         pheno5_2012=NLB033Z_2,
         sex, mob, yob, moi_2012, yoi_2012)

extraversion_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno1_2014=OLB031A,
         pheno2_2014=OLB031F,
         pheno3_2014=OLB031J,
         pheno4_2014=OLB031U,
         pheno5_2014=OLB031Z_2,
         sex, mob, yob, moi_2014, yoi_2014)

extraversion_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno1_2016=PLB031A,
         pheno2_2016=PLB031F,
         pheno3_2016=PLB031J,
         pheno4_2016=PLB031U,
         pheno5_2016=PLB031Z_2,
         sex, mob, yob, moi_2016, yoi_2016)

# transform data for residualising
extraversion <- extraversion_2006 %>%
  full_join(extraversion_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(extraversion_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(extraversion_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(extraversion_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(extraversion_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016,
         pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016,
         pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016,
         pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016,
         pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016,
         sex,
         mob,
         yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno1_", seq(2006, 2016, 2)),
         paste0("pheno2_", seq(2006, 2016, 2)),
         paste0("pheno3_", seq(2006, 2016, 2)),
         paste0("pheno4_", seq(2006, 2016, 2)),
         paste0("pheno5_", seq(2006, 2016, 2)),
         paste0("moi_",    seq(2006, 2016, 2)),
         paste0("yoi_",    seq(2006, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:5)) %>%
  arrange(HHID, PN, col, pheno_number) %>%
  group_by(HHID, PN, col) %>%
  mutate(
    pheno_count = 5 - sum(is.na(pheno_value)),
    pheno = case_when(pheno_count > 2 ~ mean(pheno_value, na.rm = T)),
    dob = yob + (mob/12),
    doi = yoi + (moi/12),
    age = doi - dob,
    age2 = age^2,
    male = 2 - sex,
    male_age = male * age,
    male_age2 = male * age2
  ) %>%
  ungroup() %>%
  group_by(HHID, PN, col, age, age2, male, male_age, male_age2) %>%
  summarise(
    pheno = 5 - mean(pheno), # reverse code
    pheno_count = mean(pheno_count)
  ) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=extraversion, average=T, age_residualise=T, name="EXTRA")


########################################################
####### Construct phenotype: family satisfaction #######
########################################################

# select and rename relevant variables
famsat_1992 <- data_1992 %>%
  select(HHID, PN,
         pheno_1992=V2615,
         sex, mob, yob, moi_1992, yoi_1992) %>%
  filter(pheno_1992 != 0)

famsat_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LLB039D,
         sex, mob, yob, moi_2008, yoi_2008)

famsat_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MLB039D,
         sex, mob, yob, moi_2010, yoi_2010)

famsat_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NLB039D,
         sex, mob, yob, moi_2012, yoi_2012)

famsat_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OLB034D,
         sex, mob, yob, moi_2014, yoi_2014)

famsat_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PLB034D,
         sex, mob, yob, moi_2016, yoi_2016)

# transform data for residualising
famsat <- famsat_1992 %>%
  full_join(famsat_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(famsat_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(famsat_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(famsat_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(famsat_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_1992, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_1992, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_1992, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1992, seq(2008, 2016, 2))),
         paste0("moi_",   c(1992, seq(2008, 2016, 2))),
         paste0("yoi_",   c(1992, seq(2008, 2016, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==5~0, pheno==4~1, pheno==3~2, pheno==2~3, pheno==1~4)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=famsat, average=T, age_residualise=T, name="FAMSAT")


########################################################
###### Construct phenotype: financial satisfaction #####
########################################################

# select and rename relevant variables
finsat_2004 <- data_2004 %>%
  select(HHID, PN,
    pheno_2004=JLB529A,
    sex, mob, yob, moi_2004, yoi_2004) %>%
    mutate(wave_2004=2004)

finsat_2006 <- data_2006 %>%
  select(HHID, PN,
    pheno_2006=KLB039A,
    sex, mob, yob, moi_2006, yoi_2006) %>%
    mutate(wave_2006=2006)

finsat_2008 <- data_2008 %>%
  select(HHID, PN,
    pheno_2008=LLB039E,
    sex, mob, yob, moi_2008, yoi_2008) %>%
    mutate(wave_2008=2008)

finsat_2010 <- data_2010 %>%
  select(HHID, PN,
    pheno_2010=MLB039E,
    sex, mob, yob, moi_2010, yoi_2010) %>%
    mutate(wave_2010=2010)

finsat_2012 <- data_2012 %>%
  select(HHID, PN,
    pheno_2012=NLB039E,
    sex, mob, yob, moi_2012, yoi_2012) %>%
    mutate(wave_2012=2012)

finsat_2014 <- data_2014 %>%
  select(HHID, PN,
    pheno_2014=OLB034E,
    sex, mob, yob, moi_2014, yoi_2014) %>%
    mutate(wave_2014=2014)

finsat_2016 <- data_2016 %>%
  select(HHID, PN,
    pheno_2016=PLB034E,
    sex, mob, yob, moi_2016, yoi_2016) %>%
    mutate(wave_2016=2016)

# transform data for residualising
finsat <- finsat_2004 %>%
  full_join(finsat_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(finsat_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(finsat_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(finsat_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(finsat_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(finsat_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2004, pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016,
         wave_2004, wave_2006, wave_2008, wave_2010, wave_2012, wave_2014, wave_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2004, 2016, 2)),
         paste0("wave_",  seq(2004, 2016, 2)),
         paste0("moi_",   seq(2004, 2016, 2)),
         paste0("yoi_",   seq(2004, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(wave>=2008 & pheno==5 ~ 0,
                           wave>=2008 & pheno==4 ~ 1,
                           wave>=2008 & pheno==3 ~ 2,
                           wave>=2008 & pheno==2 ~ 3,
                           wave>=2008 & pheno==1 ~ 4,
                           # ordering of financial satisfaction answers changes
                           wave<2008  & pheno==5 ~ 4,
                           wave<2008  & pheno==4 ~ 3,
                           wave<2008  & pheno==3 ~ 2,
                           wave<2008  & pheno==2 ~ 1,
                           wave<2008  & pheno==1 ~ 0)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=finsat, average=T, age_residualise=T, name="FINSAT")


########################################################
####### Construct phenotype: friend satisfaction #######
########################################################

# select and rename relevant variables
friendsat_1992 <- data_1992 %>%
  select(HHID, PN,
    pheno_1992=V2612,
    sex, mob, yob, moi_1992, yoi_1992) %>%
    mutate(wave_1992=1992)

# transform data for residualising
friendsat <- friendsat_1992 %>%
  select(HHID, PN,
         pheno = pheno_1992,
         sex,
         mob,
         yob,
         moi = moi_1992,
         yoi = yoi_1992,
         wave = wave_1992) %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==5~0, pheno==4~1, pheno==3~2, pheno==2~3, pheno==1~4)) %>%
  select(HHID, PN, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=friendsat, average=F, age_residualise=T, name="FRIENDSAT")

#################################################################################################
##################### Construct phenotype: Health Condition (Coronary Artery Diseas) ############
################################################################################################

# select and rename relevant variables
heartDisease_1995 <- data_1995 %>%
  select(HHID, PN,
         pheno_1995=D828,
         sex, mob, yob, moi_1995, yoi_1995)

heartDisease_1996 <- data_1996 %>%
  select(HHID, PN,
         pheno_1996=E828,
         sex, mob, yob, moi_1996, yoi_1996)

heartDisease_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=F1156,
         sex, mob, yob, moi_1998, yoi_1998)

heartDisease_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G1289,
         sex, mob, yob, moi_2000, yoi_2000)

heartDisease_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC036,
         sex, mob, yob, moi_2002, yoi_2002)

heartDisease_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC036,
         sex, mob, yob, moi_2004, yoi_2004)

heartDisease_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KC036,
         sex, mob, yob, moi_2006, yoi_2006)

heartDisease_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LC036,
         sex, mob, yob, moi_2008, yoi_2008)

heartDisease_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC036,
         sex, mob, yob, moi_2010, yoi_2010)

heartDisease_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC036,
         sex, mob, yob, moi_2012, yoi_2012)

heartDisease_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC036,
         sex, mob, yob, moi_2014, yoi_2014)

heartDisease_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC036,
         sex, mob, yob, moi_2016, yoi_2016)

heartDisease_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC036,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
heartDisease <- heartDisease_1995 %>%
  full_join(heartDisease_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(heartDisease_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,

         sex,
         mob,
         yob,

         moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,

         yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1995, seq(1996, 2018, 2))),
         paste0("moi_",   c(1995, seq(1996, 2018, 2))),
         paste0("yoi_",   c(1995, seq(1996, 2018, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  filter(!all(is.na(pheno))) %>%  # Exclude groups where all pheno values are NA
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, wave=col, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

  # Check pheno generated correclty 
   frequency_table <- table(heartDisease$pheno)
   print(frequency_table)

  # Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
  heartDisease <- heartDisease %>%
    filter(HHID %in% ancestry_hhids)  

  # Calculate summary statistics for the 'pheno' column
  summary_stats <- heartDisease %>%
    summarise(
      Min = min(pheno, na.rm = TRUE),
      Max = max(pheno, na.rm = TRUE),
      Mean = mean(pheno, na.rm = TRUE),
      Median = median(pheno, na.rm = TRUE),
      SD = sd(pheno, na.rm = TRUE),
      N = n_distinct(HHID, PN)  # Number of unique IDs
    )

  # Print the summary statistics
  print("Summary Statistics of pheno column in heartDisease:")
  print(summary_stats)

  # View the first few rows of the processed data
  print("First few rows of processed data:")
  print(head(heartDisease)) 

# residualise, average, save
residualise.average.save(data=heartDisease, average=F, age_residualise=F, name="CAD")


###########################################################################################################
##################### Construct phenotype: Health Condition (Coronary Artery Diseas - HARDCAD) ############
###########################################################################################################

# select and rename relevant variables - HARDCAD defined based on heart-attach and angina
hardcad_1992 <- data_1992 %>% 
  select(HHID, PN, 
         heartattack_1992=V407, angina_1992=V409, 
         sex, mob, yob, moi_1992, yoi_1992)
str(hardcad_1992)

hardcad_1993 <- data_1993 %>% 
  select(HHID, PN, 
         heartattack_1993=V244, angina_1993=V245, 
         sex, mob, yob, moi_1993, yoi_1993)
str(hardcad_1993)

hardcad_1994 <- data_1994 %>% 
  select(HHID, PN, 
         heartattack_1994=W369, angina_1994=W371, 
         sex, mob, yob, moi_1994, yoi_1994)
 str(hardcad_1994)        

hardcad_1995 <- data_1995 %>% 
  select(HHID, PN, 
         heartattack_1995=D834, angina_1995=D840, 
         sex, mob, yob, moi_1995, yoi_1995)
str(hardcad_1995)

hardcad_1996 <- data_1996 %>% 
  select(HHID, PN, 
         heartattack_1996=E834, angina_1996=E840, 
         sex, mob, yob, moi_1996, yoi_1996)
str(hardcad_1996)

hardcad_1998 <- data_1998 %>% 
  select(HHID, PN, 
         heartattack_1998=F1162, angina_1998=F1168, 
         sex, mob, yob, moi_1998, yoi_1998)
str(hardcad_1998)

hardcad_2000 <- data_2000 %>% 
  select(HHID, PN, 
         heartattack_2000=G1295, angina_2000=G1301, 
         sex, mob, yob, moi_2000, yoi_2000)
str(hardcad_2000)

hardcad_2002 <- data_2002 %>% 
  select(HHID, PN, 
         heartattack_2002=HC040, angina_2002=HC045, 
         sex, mob, yob, moi_2002, yoi_2002)
str(hardcad_2002)

hardcad_2004 <- data_2004 %>% 
  select(HHID, PN, 
         heartattack_2004=JC040, angina_2004=JC045, 
         sex, mob, yob, moi_2004, yoi_2004)
str(hardcad_2004)

hardcad_2006 <- data_2006 %>% 
  select(HHID, PN, 
         heartattack_2006=KC040, angina_2006=KC045, 
         sex, mob, yob, moi_2006, yoi_2006)
str(hardcad_2006)

hardcad_2008 <- data_2008 %>% 
  select(HHID, PN, 
         heartattack_2008=LC040, angina_2008=LC045, 
         sex, mob, yob, moi_2008, yoi_2008)
str(hardcad_2008)

hardcad_2010 <- data_2010 %>% 
  select(HHID, PN, 
         heartattack_2010=MC040, angina_2010=MC045, 
         sex, mob, yob, moi_2010, yoi_2010)
str(hardcad_2010)

hardcad_2012 <- data_2012 %>% 
  select(HHID, PN, 
         heartattack_2012=NC040, angina_2012=NC045, 
         sex, mob, yob, moi_2012, yoi_2012)
str(hardcad_2012)

hardcad_2014 <- data_2014 %>% 
  select(HHID, PN, 
         heartattack_2014=OC040, angina_2014=OC045, 
         sex, mob, yob, moi_2014, yoi_2014)
str(hardcad_2014)

hardcad_2016 <- data_2016 %>% 
  select(HHID, PN, 
         heartattack_2016=PC040, angina_2016=PC045, 
         sex, mob, yob, moi_2016, yoi_2016)
str(hardcad_2016)

hardcad_2018 <- data_2018 %>% 
  select(HHID, PN, 
         heartattack_2018=QC040, angina_2018=QC045, 
         sex, mob, yob, moi_2018, yoi_2018)
str(hardcad_2018)

# transform data for residualising
hardcad <- hardcad_1992 %>%
  full_join(hardcad_1993, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_1994, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(hardcad_1995, by=c("HHID", "PN", "sex", "mob", "yob")) %>%  
  full_join(hardcad_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  full_join(hardcad_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>% 
  mutate(across(starts_with("heartattack_"), ~ifelse(. == 1, 1, ifelse(. == 5, 0, NA)))) %>%
  mutate(across(starts_with("angina_"), ~ifelse(. == 1, 1, ifelse(. == 5, 0, NA)))) %>%
  select(HHID, PN,
         starts_with("heartattack_"),
         starts_with("angina_"),
         starts_with("moi_"), 
         starts_with("yoi_"),
         sex, mob, yob) %>%
  gather(key="wave", value="value",
         paste0("heartattack_", c(1992, 1993, 1994, 1995, seq(1996, 2018, 2))),
         paste0("angina_", c(1992, 1993, 1994, 1995, seq(1996, 2018, 2))), 
         paste0("moi_", c(1992, 1993, 1994, 1995, seq(1996, 2018, 2))),
         paste0("yoi_", c(1992, 1993, 1994, 1995, seq(1996, 2018, 2)))) %>%
  separate("wave", into = c("var", "year"), sep = "_") %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(
            heartattack == 1 | angina == 1 ~ 1,   # If either is 1, set to 1
            heartattack == 0 & angina == 0 ~ 0,   # If both are 0, set to 0
            TRUE ~ NA_real_                        # All other cases, set to NA
          )) %>%
  group_by(HHID, PN) %>%
  filter(!all(is.na(pheno))) %>%  # Exclude groups where all pheno values are NA
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, wave=year, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Tabulate values in the HARDCAD column
print(table_hardcad <- table(hardcad$pheno, useNA = "always"))

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
hardcad <- hardcad %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- hardcad %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in hardcad:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(hardcad)) 

# residualise, average, save
residualise.average.save(data=hardcad, average=F, age_residualise=F, name="HARDCAD")


##########################################################################
####### Construct phenotype: high cholesterol (used for BL_CHOL) ##########
###########################################################################

# select and rename relevant variables
highchol_1992 <- data_1992 %>%
  select(HHID, PN,
    pheno_1992=V435,
    sex, mob, yob, moi_1992, yoi_1992) %>%
    mutate(wave_1992=1992)
print(head(highchol_1992$hhidpn))

# transform data for residualising
highchol <- highchol_1992 %>%
  select(HHID, PN,
         pheno = pheno_1992,
         sex,
         mob,
         yob,
         moi = moi_1992,
         yoi = yoi_1992,
         wave = wave_1992) %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==5~0, pheno==1~1)) %>%
  select(HHID, PN, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Calculate summary statistics for the 'pheno' column
summary_stats <- highchol %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in highchol:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(highchol)) 

# residualise, average, save
residualise.average.save(data=highchol, average=F, age_residualise=T, name="BL_CHOL")
#residualise.average.save(data=highchol, average=F, age_residualise=T, name="BL_HDL")
#residualise.average.save(data=highchol, average=F, age_residualise=T, name="BL_LDL")
#residualise.average.save(data=highchol, average=F, age_residualise=T, name="BL_TRYG")
residualise.average.save(data=highchol, average=F, age_residualise=T, name="BL_nonHDL")

########################################################
############# Construct phenotype: Insomnia ############
########################################################

# select and rename relevant variables
insomnia_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC083,
         sex, mob, yob, moi_2002, yoi_2002)

insomnia_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC083,
         sex, mob, yob, moi_2004, yoi_2004)

insomnia_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KC083,
         sex, mob, yob, moi_2006, yoi_2006)

insomnia_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LC083,
         sex, mob, yob, moi_2008, yoi_2008)

insomnia_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC083,
         sex, mob, yob, moi_2010, yoi_2010)

insomnia_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC083,
         sex, mob, yob, moi_2012, yoi_2012)

insomnia_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC083,
         sex, mob, yob, moi_2014, yoi_2014)

insomnia_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC083,
         sex, mob, yob, moi_2016, yoi_2016)

insomnia_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC083,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
insomnia <- insomnia_2002 %>%
  full_join(insomnia_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(insomnia_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,

         sex,
         mob,
         yob,

         moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,

         yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(seq(2002, 2018, 2))),
         paste0("moi_",   c(seq(2002, 2018, 2))),
         paste0("yoi_",   c(seq(2002, 2018, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno<=2 ~ 1, pheno==3 ~ 0)) %>%
         #those trouble falling asleep most of the times and sometimes ==1; those rarely or never==0 
  group_by(HHID, PN) %>%
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, wave=col, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
insomnia <- insomnia %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- insomnia %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in insomnia:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(insomnia)) 

# residualise, average, save
residualise.average.save(data=insomnia, average=F, age_residualise=F, name="INSOMNIA")


########################################################
############# Construct phenotype: migraine ############
########################################################

# select and rename relevant variables
migraine_1995 <- data_1995 %>%
  select(HHID, PN,
         pheno_1995=D969,
         sex, mob, yob, moi_1995, yoi_1995)

migraine_1996 <- data_1996 %>%
  select(HHID, PN,
         pheno_1996=E971,
         sex, mob, yob, moi_1996, yoi_1996)

migraine_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=F1308,
         sex, mob, yob, moi_1998, yoi_1998)

migraine_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G1441,
         sex, mob, yob, moi_2000, yoi_2000)

migraine_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=HC147,
         sex, mob, yob, moi_2002, yoi_2002)

migraine_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JC147,
         sex, mob, yob, moi_2004, yoi_2004)

migraine_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2006=LC147,
         sex, mob, yob, moi_2008, yoi_2008)

migraine_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2008=KC147,
         sex, mob, yob, moi_2006, yoi_2006)

migraine_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MC147,
         sex, mob, yob, moi_2010, yoi_2010)

migraine_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NC147,
         sex, mob, yob, moi_2012, yoi_2012)

migraine_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OC147,
         sex, mob, yob, moi_2014, yoi_2014)

migraine_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PC147,
         sex, mob, yob, moi_2016, yoi_2016)

migraine_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QC147,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
migraine <- migraine_1995 %>%
  full_join(migraine_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(migraine_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,

         sex,
         mob,
         yob,

         moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,

         yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1995, seq(1996, 2018, 2))),
         paste0("moi_",   c(1995, seq(1996, 2018, 2))),
         paste0("yoi_",   c(1995, seq(1996, 2018, 2)))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = max(pheno, na.rm=T),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, wave=col, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
migraine <- migraine %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- migraine %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in migraine:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(migraine)) 

# residualise, average, save
residualise.average.save(data=migraine, average=F, age_residualise=F, name="MIGRAINE")


########################################################
######### Construct phenotype: nearsightedness #########
########################################################

# select and rename relevant variables
nearsighted_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G6893,
         sex, mob, yob, moi_2000, yoi_2000)

# transform data for residualising
nearsighted <- nearsighted_2000 %>%
  select(HHID, PN,
         pheno_2000,
         sex, mob, yob, moi_2000, yoi_2000) %>%
  gather(key="wave", value="value",
         "pheno_2000",
         "moi_2000", "yoi_2000") %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==1 ~ 0, pheno==5 ~ 1)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
nearsighted <- nearsighted %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- nearsighted %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in nearsighted:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(nearsighted)) 

# residualise, average, save
residualise.average.save(data=nearsighted, average=F, age_residualise=T, name="NEARSIGHTED")


########################################################
########### Construct phenotype: neuroticism ###########
########################################################

# select and rename relevant variables
neuro_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno1_2006=KLB033C,
         pheno2_2006=KLB033G,
         pheno3_2006=KLB033K,
         pheno4_2006=KLB033P,
         sex, mob, yob, moi_2006, yoi_2006)

neuro_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno1_2008=LLB033C,
         pheno2_2008=LLB033G,
         pheno3_2008=LLB033K,
         pheno4_2008=LLB033P,
         sex, mob, yob, moi_2008, yoi_2008)

neuro_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno1_2010=MLB033D,
         pheno2_2010=MLB033H,
         pheno3_2010=MLB033L,
         pheno4_2010=MLB033Q,
         sex, mob, yob, moi_2010, yoi_2010)

neuro_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno1_2012=NLB033D,
         pheno2_2012=NLB033H,
         pheno3_2012=NLB033L,
         pheno4_2012=NLB033Q,
         sex, mob, yob, moi_2012, yoi_2012)

neuro_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno1_2014=OLB031D,
         pheno2_2014=OLB031H,
         pheno3_2014=OLB031L,
         pheno4_2014=OLB031Q,
         sex, mob, yob, moi_2014, yoi_2014)

neuro_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno1_2016=PLB031D,
         pheno2_2016=PLB031H,
         pheno3_2016=PLB031L,
         pheno4_2016=PLB031Q,
         sex, mob, yob, moi_2016, yoi_2016)

neuro_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno1_2018=QLB031D,
         pheno2_2018=QLB031H,
         pheno3_2018=QLB031L,
         pheno4_2018=QLB031Q,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
neuro <- neuro_2006 %>%
  full_join(neuro_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno1_2006, pheno1_2008, pheno1_2010,
         pheno1_2012, pheno1_2014, pheno1_2016, pheno1_2018, 

         pheno2_2006, pheno2_2008, pheno2_2010,
         pheno2_2012, pheno2_2014, pheno2_2016, pheno2_2018,

         pheno3_2006, pheno3_2008, pheno3_2010,
         pheno3_2012, pheno3_2014, pheno3_2016, pheno3_2018,  

         pheno4_2006, pheno4_2008, pheno4_2010,
         pheno4_2012, pheno4_2014, pheno4_2016, pheno4_2018, 

         sex,
         mob,
         yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno1_", seq(2006, 2018, 2)),
         paste0("pheno2_", seq(2006, 2018, 2)),
         paste0("pheno3_", seq(2006, 2018, 2)),
         paste0("pheno4_", seq(2006, 2018, 2)),
         paste0("moi_",    seq(2006, 2018, 2)),
         paste0("yoi_",    seq(2006, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(
    pheno1 = 5 - pheno1,
    pheno2 = 5 - pheno2,
    pheno3 = 5 - pheno3
  ) %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:4)) %>%
  arrange(HHID, PN, col, pheno_number) %>%
  group_by(HHID, PN, col) %>%
  mutate(
    pheno_count = 4 - sum(is.na(pheno_value)),
    pheno = case_when(pheno_count > 1 ~ mean(pheno_value, na.rm = T)),
    dob = yob + (mob/12),
    doi = yoi + (moi/12),
    age = doi - dob,
    age2 = age^2,
    male = 2 - sex,
    male_age = male * age,
    male_age2 = male * age2
  ) %>%
  ungroup() %>%
  group_by(HHID, PN, col, age, age2, male, male_age, male_age2) %>%
  summarise(
    pheno = 5 - mean(pheno), # reverse code
    pheno_count = mean(pheno_count)
  ) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=neuro, average=T, age_residualise=T, name="NEURO")



########################################################
############# Construct phenotype: openness ############
########################################################

# select and rename relevant variables
open_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno1_2006=KLB033L,
         pheno2_2006=KLB033N,
         pheno3_2006=KLB033Q,
         pheno4_2006=KLB033R,
         pheno5_2006=KLB033U,
         pheno6_2006=KLB033X,
         pheno7_2006=KLB033Y,
         sex, mob, yob, moi_2006, yoi_2006)

open_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno1_2008=LLB033L,
         pheno2_2008=LLB033N,
         pheno3_2008=LLB033Q,
         pheno4_2008=LLB033R,
         pheno5_2008=LLB033U,
         pheno6_2008=LLB033X,
         pheno7_2008=LLB033Y,
         sex, mob, yob, moi_2008, yoi_2008)

open_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno1_2010=MLB033M,
         pheno2_2010=MLB033O,
         pheno3_2010=MLB033S,
         pheno4_2010=MLB033T,
         pheno5_2010=MLB033W,
         pheno6_2010=MLB033Z_3,
         pheno7_2010=MLB033Z_4,
         sex, mob, yob, moi_2010, yoi_2010)

open_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno1_2012=NLB033M,
         pheno2_2012=NLB033O,
         pheno3_2012=NLB033S,
         pheno4_2012=NLB033T,
         pheno5_2012=NLB033W,
         pheno6_2012=NLB033Z_3,
         pheno7_2012=NLB033Z_4,
         sex, mob, yob, moi_2012, yoi_2012)

open_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno1_2014=OLB031M,
         pheno2_2014=OLB031O,
         pheno3_2014=OLB031S,
         pheno4_2014=OLB031T,
         pheno5_2014=OLB031W,
         pheno6_2014=OLB031Z_3,
         pheno7_2014=OLB031Z_4,
         sex, mob, yob, moi_2014, yoi_2014)

open_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno1_2016=PLB031M,
         pheno2_2016=PLB031O,
         pheno3_2016=PLB031S,
         pheno4_2016=PLB031T,
         pheno5_2016=PLB031W,
         pheno6_2016=PLB031Z_3,
         pheno7_2016=PLB031Z_4,
         sex, mob, yob, moi_2016, yoi_2016)

open_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno1_2018=QLB031M,
         pheno2_2018=QLB031O,
         pheno3_2018=QLB031S,
         pheno4_2018=QLB031T,
         pheno5_2018=QLB031W,
         pheno6_2018=QLB031Z_3,
         pheno7_2018=QLB031Z_4,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
open <- open_2006 %>%
  full_join(open_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016, pheno1_2018,
         pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016, pheno2_2018,
         pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016, pheno3_2018,
         pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016, pheno4_2018,
         pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016, pheno5_2018,
         pheno6_2006, pheno6_2008, pheno6_2010, pheno6_2012, pheno6_2014, pheno6_2016, pheno6_2018,
         pheno7_2006, pheno7_2008, pheno7_2010, pheno7_2012, pheno7_2014, pheno7_2016, pheno7_2018,
         sex,
         mob,
         yob,
         moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno1_", seq(2006, 2018, 2)),
         paste0("pheno2_", seq(2006, 2018, 2)),
         paste0("pheno3_", seq(2006, 2018, 2)),
         paste0("pheno4_", seq(2006, 2018, 2)),
         paste0("pheno5_", seq(2006, 2018, 2)),
         paste0("pheno6_", seq(2006, 2018, 2)),
         paste0("pheno7_", seq(2006, 2018, 2)),
         paste0("moi_",    seq(2006, 2018, 2)),
         paste0("yoi_",    seq(2006, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  gather(key = "pheno_number", value = "pheno_value", paste0("pheno", 1:7)) %>%
  arrange(HHID, PN, col, pheno_number) %>%
  group_by(HHID, PN, col) %>%
  mutate(
    pheno_count = 7 - sum(is.na(pheno_value)),
    pheno = case_when(pheno_count > 3 ~ mean(pheno_value, na.rm = T)),
    dob = yob + (mob/12),
    doi = yoi + (moi/12),
    age = doi - dob,
    age2 = age^2,
    male = 2 - sex,
    male_age = male * age,
    male_age2 = male * age2
  ) %>%
  ungroup() %>%
  group_by(HHID, PN, col, age, age2, male, male_age, male_age2) %>%
  summarise(
    pheno = 5 - mean(pheno), # reverse code
    pheno_count = mean(pheno_count)
  ) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=open, average=T, age_residualise=T, name="OPEN")


########################################################
########### Construct phenotype: religiosity ###########
########################################################

# select and rename relevant variables

# About how often have you attended religious services during the past year?
#1.      More than once a week
#2.      Once a week
#3.      Two or three times a month
#4.      One or more times a year
#5.      Not at all
# 7 other
#8.      DK
#9.      NA
#0.      Inap, 90, 98-99 in V214
religiosity_1992 <- data_1992 %>%
  select(HHID, PN,
         pheno_1992=V215,
         sex, mob, yob, moi_1992, yoi_1992) %>%
  mutate(wave=1992)

religiosity_1994 <- data_1994 %>%
  select(HHID, PN,
         pheno_1994=W227,
         sex, mob, yob, moi_1994, yoi_1994) %>%
  mutate(wave=1994)

religiosity_1995 <- data_1995 %>%
  select(HHID, PN,
         pheno_1995=D736,
         sex, mob, yob, moi_1995, yoi_1995) %>%
  mutate(wave=1995)

religiosity_1996 <- data_1996 %>%
  select(HHID, PN,
         pheno_1996=P754,
         sex, mob, yob, moi_1996, yoi_1996) %>%
  mutate(wave=1996) # from 1996 EXIT survey

religiosity_1998 <- data_1998 %>%
  select(HHID, PN,
         pheno_1998=Q1056,
         sex, mob, yob, moi_1998, yoi_1998) %>%
  mutate(wave_1998=1998) # from 1998 EXIT survey
#1. DAILY
#2. AT LEAST 2 OR 3 TIMES A WEEK
#3. AT LEAST ONCE A WEEK
#4. AT LEAST ONCE A MONTH
#5. LESS THAN ONCE A MONTH
#8. DK (don't know); NA (not ascertained)
#9. RF (refused)
religiosity_2000 <- data_2000 %>%
  select(HHID, PN,
         pheno_2000=G6773,
         sex, mob, yob, moi_2000, yoi_2000) %>%
  mutate(wave_2000=2000)
#          2217           1.  MORE THAN ONCE A WEEK
#          4124           2.  ONCE A WEEK
#          2105           3.  TWO OR THREE TIMES A MONTH
#          3253           4.  ONE OR MORE TIMES A YEAR
#          4465           5.  NOT AT ALL
#           117           8.  DK (Dont Know); NA (Not Ascertained)
#            15           9.  RF (Refused)
#                     Blank.  INAP (Inapplicable); Partial Interview
religiosity_2002 <- data_2002 %>%
  select(HHID, PN,
         pheno_2002=SB082,
         sex, mob, yob, moi_2002, yoi_2002) %>%
  mutate(wave_2002=2002) # from 2002 EXIT survey

religiosity_2004 <- data_2004 %>%
  select(HHID, PN,
         pheno_2004=JB082,
         sex, mob, yob, moi_2004, yoi_2004) %>%
  mutate(wave_2004=2004)

religiosity_2006 <- data_2006 %>%
  select(HHID, PN,
         pheno_2006=KB082,
         sex, mob, yob, moi_2006, yoi_2006) %>%
  mutate(wave_2006=2006)

religiosity_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LB082,
         sex, mob, yob, moi_2008, yoi_2008) %>%
  mutate(wave_2008=2008)

religiosity_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MB082,
         sex, mob, yob, moi_2010, yoi_2010) %>%
  mutate(wave_2010=2010)

religiosity_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NB082,
         sex, mob, yob, moi_2012, yoi_2012) %>%
  mutate(wave_2012=2012)

religiosity_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OB082,
         sex, mob, yob, moi_2014, yoi_2014) %>%
  mutate(wave_2014=2014)

religiosity_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PB082,
         sex, mob, yob, moi_2016, yoi_2016) %>%
  mutate(wave_2016=2016)

religiosity_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QB082,
         sex, mob, yob, moi_2018, yoi_2018) %>%
  mutate(wave_2018=2018)

### coding
# 0 not at all
# 1 one or more times a year OR less than once a month (inapplicable left blank)
# 2 two or three times a month OR at least once a month
# 3 once a week
# 4 two or three times a week OR more than once a week
# 5 daily

# transform data for residualising
religiosity <- religiosity_1992 %>%
  full_join(religiosity_1994, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_1995, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#  full_join(religiosity_1996, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#  full_join(religiosity_1998, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2000, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
#  full_join(religiosity_2002, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(religiosity_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2004, pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,
         wave_2004, wave_2006, wave_2008, wave_2010, wave_2012, wave_2014, wave_2016, wave_2018,
         sex,
         mob,
         yob,
         moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2004, 2018, 2)),
         paste0("wave_",  seq(2004, 2018, 2)),
         paste0("moi_",   seq(2004, 2018, 2)),
         paste0("yoi_",   seq(2004, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(wave>2000  & pheno==5 ~ 0,
                           wave>2000  & pheno==4 ~ 1,
                           wave>2000  & pheno==3 ~ 2,
                           wave>2000  & pheno==2 ~ 3,
                           wave>2000  & pheno==1 ~ 4,

                           wave==2000 & pheno==5 ~ 1,
                           wave==2000 & pheno==4 ~ 2,
                           wave==2000 & pheno==3 ~ 3,
                           wave==2000 & pheno==2 ~ 4,
                           wave==2000 & pheno==1 ~ 5,

                           wave<2000  & pheno==5 ~ 0,
                           wave<2000  & pheno==4 ~ 1,
                           wave<2000  & pheno==3 ~ 2,
                           wave<2000  & pheno==2 ~ 3,
                           wave<2000  & pheno==1 ~ 4)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=religiosity, average=T, age_residualise=T, name="RELIGATT")

#############################################################
######### Construct phenotype: schizophrenia ################
#############################################################

# select and rename relevant variables
schizophrenia_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OV354,
         sex, mob, yob, moi_2014, yoi_2014)

# transform data for residualising
schizophrenia <- schizophrenia_2014 %>%
  select(HHID, PN, pheno_2014, sex, mob, yob, moi_2014, yoi_2014) %>%
  gather(key="wave", value="value",
         paste0("pheno_", 2014),
         paste0("moi_",   2014),
         paste0("yoi_",   2014)) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         dob2 = dob^2,
         male = 2 - sex,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = case_when(pheno==1 ~ 1, pheno==5 ~ 0)) %>%
  group_by(HHID, PN) %>%
  mutate(pheno = ifelse(!is.na(pheno) & pheno < 8, max(pheno, na.rm = TRUE), NA),
         rn = row_number()) %>%
  ungroup() %>%
  filter(rn==1 & pheno!="-Inf") %>%
  select(HHID, PN, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# Filter the breastCancer data frame for individuals with HHID in the requested ancestry HHIDs list
schizophrenia <- schizophrenia %>%
  filter(HHID %in% ancestry_hhids)  

# Calculate summary statistics for the 'pheno' column
summary_stats <- schizophrenia %>%
  summarise(
    Min = min(pheno, na.rm = TRUE),
    Max = max(pheno, na.rm = TRUE),
    Mean = mean(pheno, na.rm = TRUE),
    Median = median(pheno, na.rm = TRUE),
    SD = sd(pheno, na.rm = TRUE),
    N = n_distinct(HHID, PN)  # Number of unique IDs
  )

# Print the summary statistics
print("Summary Statistics of pheno column in schizophrenia:")
print(summary_stats)

# View the first few rows of the processed data
print("First few rows of processed data:")
print(head(schizophrenia)) 

# residualise, average, save
residualise.average.save(data=schizophrenia, average=F, age_residualise=F, name="SCZ")


########################################################
######## Construct phenotype: self-rated health ########
########################################################

# select and rename relevant variables
selfhealth_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LLB039F,
         sex, mob, yob, moi_2008, yoi_2008)

selfhealth_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MLB039G,
         sex, mob, yob, moi_2010, yoi_2010)

selfhealth_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NLB039G,
         sex, mob, yob, moi_2012, yoi_2012)

selfhealth_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OLB034G,
         sex, mob, yob, moi_2014, yoi_2014)

selfhealth_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PLB034G,
         sex, mob, yob, moi_2016, yoi_2016)

selfhealth_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QLB034G,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
selfhealth <- selfhealth_2008 %>%
  full_join(selfhealth_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016, pheno_2018,
         sex,
         mob,
         yob,
         moi_2008, moi_2010, moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2008, 2018, 2)),
         paste0("moi_",   seq(2008, 2018, 2)),
         paste0("yoi_",   seq(2008, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==5~0, pheno==4~1, pheno==3~2, pheno==2~3, pheno==1~4)) %>%
  select(HHID, PN, wave=col, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=selfhealth, average=T, age_residualise=T, name="SELFHEALTH")


########################################################
######## Construct phenotype: work satisfaction ########
########################################################

# select and rename relevant variables
worksat_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NLB084A,
         sex, mob, yob, moi_2012, yoi_2012)

worksat_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OLB076,
         sex, mob, yob, moi_2014, yoi_2014)

worksat_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PLB076,
         sex, mob, yob, moi_2016, yoi_2016)

worksat_2018 <- data_2018 %>%
  select(HHID, PN,
         pheno_2018=QLB076,
         sex, mob, yob, moi_2018, yoi_2018)

# transform data for residualising
worksat <- worksat_2012 %>%
  full_join(worksat_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(worksat_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(worksat_2018, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2012, pheno_2014, pheno_2016, pheno_2018,
         sex,
         mob,
         yob,
         moi_2012, moi_2014, moi_2016, moi_2018,
         yoi_2012, yoi_2014, yoi_2016, yoi_2018) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2012, 2018, 2)),
         paste0("moi_",   seq(2012, 2018, 2)),
         paste0("yoi_",   seq(2012, 2018, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%

  mutate(dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno==1~0, pheno==2~1, pheno==3~2, pheno==4~3)) %>%
  select(HHID, PN, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=worksat, average=F, age_residualise=T, name="WORKSAT")
