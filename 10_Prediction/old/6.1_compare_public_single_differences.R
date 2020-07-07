#----------------------------------------------------------------------------------#
# Compares differences between public and single-trait scores
# Date: 04/21/2020
# Author: Joel Becker

# Notes:
#  TODO: not sure how NEBPOOLED supposed to be working
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# source and set directory
source("directory_paths.R")
setwd(maindir)


########################################################
####################### Load data ######################
########################################################

results_wd <- paste0(getwd(), "/jbecker/output/prediction/")

bootstraps_HRS_single_path <- paste0(results_wd, "HRS_bootstraps_single.txt")
bootstraps_WLS_single_path <- paste0(results_wd, "WLS_bootstraps_single.txt")
bootstraps_HRS_public_path <- paste0(results_wd, "HRS_bootstraps_public.txt")
bootstraps_WLS_public_path <- paste0(results_wd, "WLS_bootstraps_public.txt")

bootstraps_HRS_single_data <- fread(bootstraps_HRS_single_path)
bootstraps_WLS_single_data <- fread(bootstraps_WLS_single_path)
bootstraps_HRS_public_data <- fread(bootstraps_HRS_public_path)

## remove HRS public NEBwomen results, which have overlap
#bootstraps_HRS_public_NEBwomen <- filter(bootstraps_HRS_public_data, phenotype == "NEBwomen")
#bootstraps_HRS_public_NEBwomen$inc_r2_public <- NA
#bootstraps_HRS_public_data <- bootstraps_HRS_public_data %>%
#  filter(phenotype != "NEBwomen") %>%
#  rbind(., bootstraps_HRS_public_NEBwomen)
#
## remove HRS public EA results, which have overlap
#bootstraps_HRS_public_EA <- filter(bootstraps_HRS_public_data, phenotype == "EA")
#bootstraps_HRS_public_EA$inc_r2_public <- NA
#bootstraps_HRS_public_data <- bootstraps_HRS_public_data %>%
#  filter(phenotype != "EA") %>%
#  rbind(., bootstraps_HRS_public_EA)
#
## remove HRS public EVERSMOKE results, which have overlap
#bootstraps_HRS_public_EVERSMOKE <- filter(bootstraps_HRS_public_data, phenotype == "EVERSMOKE")
#bootstraps_HRS_public_EVERSMOKE$inc_r2_public <- NA
#bootstraps_HRS_public_data <- bootstraps_HRS_public_data %>%
#  filter(phenotype != "EVERSMOKE") %>%
#  rbind(., bootstraps_HRS_public_EVERSMOKE)

bootstraps_WLS_public_data <- fread(bootstraps_WLS_public_path)

## remove WLS public EA results, which have overlap
#bootstraps_WLS_public_EA <- filter(bootstraps_WLS_public_data, phenotype == "EA")
#bootstraps_WLS_public_EA$inc_r2_public <- NA
#bootstraps_WLS_public_data <- bootstraps_WLS_public_data %>%
#  filter(phenotype != "EA") %>%
#  rbind(., bootstraps_WLS_public_EA)


########################################################
################ Rename results to match ###############
########################################################

#bootstraps_HRS_single_data <- bootstraps_HRS_single_data %>%
#  mutate(
#    phenotype = case_when(
#      phenotype == "AsthmaEczRhi" ~ "ASTECZRHI",
#      phenotype == "RELIGIOUS" ~ "RELIGIOSITY",
#      TRUE ~ phenotype
#    )
#  )
#
#bootstraps_WLS_single_data <- bootstraps_WLS_single_data %>%
#  mutate(
#    phenotype = case_when(
#      phenotype == "AsthmaEczRhi" ~ "ASTECZRHI",
#      phenotype == "AGE_FIRST_MENSES" ~ "AGEFIRSTMENSES",
#      #phenotype == "AsthmaEczRhi" ~ "NEBPOOLED",
#      phenotype == "RELIGIOUS" ~ "RELIGIOSITY",
#      TRUE ~ phenotype
#    )
#  )

bootstraps_HRS_public_data <- bootstraps_HRS_public_data %>%
  mutate(
    phenotype = case_when(
      phenotype == "NEBWOMEN" ~ "NEBwomen",
      TRUE ~ phenotype
    )
  )

bootstraps_WLS_public_data <- bootstraps_WLS_public_data %>%
  mutate(
    phenotype = case_when(
      phenotype == "NEBWOMEN" ~ "NEBwomen",
      TRUE ~ phenotype
    )
  )

# append NEBWOM single-scores data to NEBPOOLED public-scores data
#WLS_NEBPOOLED <- bootstraps_WLS_single_data %>%
#  filter(phenotype == "NEBWOM") %>%
#  mutate(phenotype = "NEBPOOLED")
#
#bootstraps_WLS_single_data <- rbind(bootstraps_WLS_single_data, WLS_NEBPOOLED)


########################################################
#################### Wrangle results ###################
########################################################

bootstraps_HRS_all_data <- bootstraps_HRS_single_data %>%
  full_join(
    bootstraps_HRS_public_data,
    by = c("iteration", "phenotype")
  ) %>%
  mutate(
    inc_r2_public = case_when(
      is.na(inc_r2_public) ~ 0,
      TRUE ~ inc_r2_public
    ),
    delta = inc_r2_single - inc_r2_public,
    cohort = "HRS"
  )

bootstraps_WLS_all_data <- bootstraps_WLS_single_data %>%
  full_join(
    bootstraps_WLS_public_data,
    by = c("iteration", "phenotype")
  ) %>%
  mutate(
    inc_r2_public = case_when(
      is.na(inc_r2_public) ~ 0,
      TRUE ~ inc_r2_public
    ),
    delta = inc_r2_single - inc_r2_public,
    cohort = "WLS"
  ) #%>%
  #filter(phenotype != "NEBWOM_UKBonly")

# select only desired phenotypes
unused_phenotypes <- c(
  "NEARSIGHTED"
)

tests_all_data <- rbind(
    bootstraps_HRS_all_data,
    bootstraps_WLS_all_data
  ) %>%
  filter(!(phenotype %in% unused_phenotypes)) %>%
  mutate(single_data_only = case_when(inc_r2_public == 0 ~ 1, TRUE ~ 0)) %>%
  group_by(phenotype, cohort) %>%
  summarise(
    single_data_only = mean(single_data_only),
    delta_mean = mean(delta),
    delta_sd = sd(delta),
    delta_lower = delta_mean - (1.96 * delta_sd),
    delta_upper = delta_mean + (1.96 * delta_sd)
  ) %>%
  # select only desired phenotypes
  filter(!(phenotype %in% unused_phenotypes)) %>%
  # remove phenotypes where cohorts overlap
  #filter(!(phenotype == "EA"        & cohort == "HRS")) %>%
  #filter(!(phenotype == "EA"        & cohort == "WLS")) %>%
  #filter(!(phenotype == "EVERSMOKE" & cohort == "HRS")) %>%
  #filter(!(phenotype == "NEBwomen"  & cohort == "HRS")) %>%
  #filter(!(phenotype == "SWB" & cohort == "HRS")) %>%
  # format digits
  mutate(
    delta_mean  = formatC(delta_mean,  digits = 2, format = "f"),
    delta_lower = formatC(delta_lower, digits = 2, format = "f"),
    delta_upper = formatC(delta_upper, digits = 2, format = "f")
  ) %>%
  select(-delta_sd) %>%
  as.data.frame(.)

tests_all_data


########################################################
###################### Save data #######################
########################################################

fwrite(tests_all_data,
  paste0(getwd(),
    "/jbecker/output/prediction/public_single_tests_r2_",
    Sys.Date(),
    ".txt")
  )

# overwrite main file
fwrite(tests_all_data,
  paste0(getwd(),
    "/jbecker/output/prediction/public_single_tests_r2.txt")
  )
