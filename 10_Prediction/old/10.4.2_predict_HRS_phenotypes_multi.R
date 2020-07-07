#----------------------------------------------------------------------------------#
# Predicts HRS phenotypes using multitrait scores
# Date: 03/03/2020
# Author: Joel Becker

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "stringr")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# source and set directory
setwd("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction")


########################################################
###################### Data paths ######################
########################################################

# scores
score_wd <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/9_Scores/multi/scores/"
score_files <- list.files(score_wd)
HRS_files_all <- score_files[grep("HRS", score_files)]
score_names <- gsub("PGS_HRS2_", "", gsub("-multi_LDpred_p1.txt", "", HRS_files_all))

# phenos
pheno_wd <- "input/HRS"
pheno_files <- list.files(pheno_wd)
pheno_names <- gsub(".pheno", "", pheno_files)

# overlapping names
both_names <- pheno_names[pheno_names %in% score_names]
score_names[!(score_names %in% both_names)]; pheno_names; both_names
pheno_names <- both_names

# scores-phenos crosswalk
score_pheno_crosswalk_path <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS/HRS_GENOTYPEV2_XREF.dta"
score_pheno_crosswalk_data <- read.dta(score_pheno_crosswalk_path) %>%
  mutate(IID = as.numeric(LOCAL_ID),
         HHID = as.numeric(HHID),
         PN = as.numeric(PN))

# PCs and related crosswalk
#PCs_crosswalk_path <- "jbecker/data/Crossref/crossref_version12.dta"
#PCs_crosswalk_data <- read.dta(PCs_crosswalk_path) %>%
#  mutate(iid = as.numeric(LOCAL_ID))

PCs_path <- "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/8_PCs/HRS2/HRS2_PCs.eigenvec"
PCs_oldnames <- paste0("V", 3:22)
PCs_newnames <- paste0("pc", 1:20)
PCs_data <- fread(PCs_path) %>%
  rename(IID = V2) %>%
  rename_at(vars(PCs_oldnames), ~ PCs_newnames)

### empty data frame with r^2
df <- data.frame(phenotype = character())

# options for prediction loop
iterations <- 1000
binary_phenotypes <- c(
  "ADHD",
  "ASTECZRHI",
  "ASTHMA",
  "EVERSMOKE",
  "HAYFEVER",
  "MIGRAINE",
  "NEARSIGHTED"
)

use_nagelkerke <- FALSE

for (i in 1:length(pheno_names)){

  local_df <- data.frame(phenotype = character(1))

  pheno <- pheno_names[i]
  local_df$phenotype <- pheno

  pheno_path <- paste0(getwd(), "/input/HRS/", pheno, ".pheno")
  score_path <- paste0(score_wd, HRS_files_all[grep(paste0(pheno, "-"), HRS_files_all, fixed=TRUE)])

  # get rid of bad paths
  score_path <- score_path[!grepl("\\*", score_path)]

  pheno_data <- fread(pheno_path)
  score_data <- fread(score_path) %>%
    rename(score = contains(substr(pheno, 1, 3)))

  if ("hhidpn" %in% colnames(pheno_data)){ # if needing to join on hhidpn variables
    score_pheno_crosswalk_data <- score_pheno_crosswalk_data %>% mutate(hhidpn = (1000 * HHID) + PN)
    data <- inner_join(score_data, score_pheno_crosswalk_data, by="IID") %>%
      inner_join(pheno_data, by="hhidpn") %>%
      inner_join(PCs_data, by="IID") %>%
      drop_na(phenotype, score)
  } else { # else, if joining on hhid and pn
    if ("hhid" %in% colnames(pheno_data)) {
      pheno_data <- pheno_data %>% mutate(HHID = hhid, PN = pn)
    }
    data <- inner_join(score_data, score_pheno_crosswalk_data, by="IID") %>%
      inner_join(pheno_data, by=c("HHID", "PN")) %>%
      #inner_join(PCs_data, by=c("HHID", "PN")) %>%
      inner_join(PCs_data, by = c("IID")) %>%
      drop_na(phenotype, score)
  }

  # get PCs
  PCs <- str_c(paste0("pc", 1:10), collapse=" + ")

  # regress with and without scores, extract r^2
  formula_no_score <- as.formula(paste0("phenotype ~ ", PCs))
  formula_with_score <- as.formula(paste0("phenotype ~ score + ", PCs))
  reg_no_score <- lm(formula_no_score, data)
  reg_with_score <- lm(formula_with_score, data)
  if (pheno %in% binary_phenotypes & use_nagelkerke) {
    # if outcome binary, calculate nagelkerke r^2

    r2_no_score <- nagelkerke(resampled_no_score_reg_fit)
    r2_with_score <- nagelkerke(resampled_with_score_reg_fit)

    local_df$r2_no_score <- extract.nagelkerke.r2(r2_no_score)
    local_df$r2_with_score <- extract.nagelkerke.r2(r2_with_score)
  } else {
    # if outcome not binary, calculate standard r^2

    local_df$r2_no_score <- summary(reg_no_score)$r.squared * 100
    local_df$r2_with_score <- summary(reg_with_score)$r.squared * 100
  }
  local_df$r2_inc <- local_df$r2_with_score - local_df$r2_no_score

  # bootstrap
  r2_inc_list <- c()

  # iterate for bootstrapping
  for (j in 1:iterations) {

    # print iteration and phenotype
    if (j %% 100 == 0){
      print(paste0("Iteration: ", j, ". Phenotype: ", pheno, "."))
    }

    # resample, regress, add r2
    resampled_data <- sample_n(data, nrow(data), replace=TRUE)
    resampled_no_score_reg_fit <- lm(formula_no_score, resampled_data)
    resampled_with_score_reg_fit <- lm(formula_with_score, resampled_data)

    if (pheno %in% binary_phenotypes & use_nagelkerke) {
      # if outcome binary, calculate nagelkerke r^2

      if (j == 1){
        print(paste0("Using Nagelkerke R^2 for phenotype: ", pheno, "."))
      }

      r2_no_score <- nagelkerke(resampled_no_score_reg_fit)
      r2_with_score <- nagelkerke(resampled_with_score_reg_fit)

      r2_no_score <- extract.nagelkerke.r2(r2_no_score)
      r2_with_score <- extract.nagelkerke.r2(r2_with_score)
    } else {
      # if outcome not binary, calculate standard r^2
      resampled_no_score_reg <- summary(resampled_no_score_reg_fit)
      resampled_with_score_reg <- summary(resampled_with_score_reg_fit)

      r2_no_score <- resampled_no_score_reg$r.squared * 100
      r2_with_score <- resampled_with_score_reg$r.squared * 100
    }
    r2_inc_list[j] <- r2_with_score - r2_no_score

  }

  # generate quantities of interest
  r2_inc_low  <- quantile(r2_inc_list, .025)[[1]]
  r2_inc_high <- quantile(r2_inc_list, .975)[[1]]
  r2_inc_mean <- mean(r2_inc_list)

  # print r2 confidence interval
  print(paste0("Incremental R^2 (95% CI): ", r2_inc_mean, " [", r2_inc_low, " : ", r2_inc_high, "]"))

  ### bootstrap function
  local_df$r2_inc_lower <- r2_inc_low
  local_df$r2_inc_upper <- r2_inc_high

  # print and record merged prediction sample N
  local_df$N <- nrow(data)
  print(paste0("N: ", nrow(data)))

  # bind local df to master df
  df <- rbind(local_df, df)

}


# overwrite main file
fwrite(
  df,
  paste0(
    getwd(),
    "/output/HRS_phenotypes_r2_multi.txt"
  )
)
