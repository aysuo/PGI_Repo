#----------------------------------------------------------------------------------#
# Finishes constructing CIDI phenotype from HRS data
# Date: 04/09/2018
# Author: Joel Becker

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "Rmpfr", "sjmisc")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# set directory
setwd("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction/")
inputDataDir = "/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS"

########################################################
####################### Load data ######################
########################################################
CIDI_path <- "tmp/DEPRESS.dta"
RAND_path <- paste0(inputDataDir,"randhrs1992_2014v2.dta")

CIDI_data <- read.dta(CIDI_path) %>%
  rename(phenotype = Z_DEP_std) %>%
  select(-contains("pc"))

fwrite(CIDI_data, paste0("output/HRS/phenotypes/", "DEP", ".pheno"))

