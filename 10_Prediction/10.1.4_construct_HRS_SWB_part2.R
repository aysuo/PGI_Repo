#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------#
# Finishes constructing SWB phenotype from HRS data
# Date: 04/16/2018
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

# set wd
setwd("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction")


########################################################
####################### Load data ######################
########################################################

SWB_path <- "tmp/SWB_all.dta"
SWB_data <- read.dta(SWB_path)


########################################################
################ Residualising function ################
########################################################

residualise <- function(data) {

  # residualise within-wave, record missings
  reg <- summary(lm(pheno ~ age + age2 + male + male_age + male_age2, data))
  sel <- which(!is.na(data$pheno)) # which observations do we have proxy phenotype for?
  data$resid <- NA
  data$resid[sel] <- reg$resid
  data$missing <- 1
  data$missing[sel] <- 0

  # standardise residuals
  mean <- mean(data$resid); sd <- sd(data$resid)
  data$std_resid <- (data$resid - mean) / sd

  return(data)

}

residualise.average.save <- function(data, name) {

  # get data column names for new df
  df <- data[FALSE,]

  # list waves, to run regressions separately
  waves <- sort(unique(data$wave))

  for (i in waves) {
    # residualise within-wave
    data_wave <- filter(data, wave==i)
    data_wave <- residualise(data=data_wave)

    # record this wave
    df <- rbind(df, data_wave)
  }

  # average residuals within-wave
  df <- df %>%
    group_by(hhid, pn) %>%
    summarise(phenotype = mean(std_resid, na.rm=TRUE)) %>%
    ungroup() %>%
    filter(!is.na(phenotype))

  # save data
  fwrite(df, paste0("tmp/",name,".pheno"))

}


########################################################
##################### Wrangle data #####################
########################################################

PA <- SWB_data %>%
  select(hhid, pn, contains("PA"), contains("sex"), contains("mob"),
         contains("yob"), contains("moi"), contains("yoi")) %>%
  gather(key="wave", value="value",
         paste0("PA_", seq(2006, 2016, 2)), paste0("sex_", seq(2006, 2016, 2)),
         paste0("mob_",  seq(2006, 2016, 2)), paste0("yob_", seq(2006, 2016, 2)),
         paste0("moi_",  seq(2006, 2016, 2)), paste0("yoi_", seq(2006, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  filter(!is.na(PA)) %>%
  filter(mob %in% 1:12 & moi %in% 1:12) %>%
  mutate(sex = as.numeric(sex),
         mob = as.numeric(mob), yob = as.numeric(yob),
         moi = as.numeric(moi), yoi = as.numeric(yoi),
         dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>% # could clean by excluding people >3 years older between interviews
  select(hhid, pn, wave=col, pheno=PA, age, age2, male, male_age, male_age2) %>%
  drop_na()

LS <- SWB_data %>%
  select(hhid, pn, contains("LS"), contains("sex"), contains("mob"),
         contains("yob"), contains("moi"), contains("yoi")) %>%
  gather(key="wave", value="value",
         paste0("LS_", seq(2006, 2016, 2)), paste0("sex_", seq(2006, 2016, 2)),
         paste0("mob_",  seq(2006, 2016, 2)), paste0("yob_", seq(2006, 2016, 2)),
         paste0("moi_",  seq(2006, 2016, 2)), paste0("yoi_", seq(2006, 2016, 2))) %>%
  separate("wave", c("var", "col")) %>%
  spread("var", "value") %>%
  filter(!is.na(LS)) %>%
  filter(mob %in% 1:12 & moi %in% 1:12) %>%
  mutate(sex = as.numeric(sex),
         mob = as.numeric(mob), yob = as.numeric(yob),
         moi = as.numeric(moi), yoi = as.numeric(yoi),
         dob = yob + (mob/12),
         doi = yoi + (moi/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2) %>% # could clean by excluding people >3 years older between interviews
  select(hhid, pn, wave=col, pheno=LS, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=PA, name="SWBPA")
residualise.average.save(data=LS, name="SWBLS")

# resave SWB as average of residualised PA and LS
SWBPA <- fread("tmp/SWBPA.pheno")
SWBLS <- fread("tmp/SWBLS.pheno")

SWB <- SWBPA %>%
  rename(PA = phenotype) %>%
  full_join(SWBLS, by = c("hhid", "pn")) %>%
  rename(LS = phenotype) %>%
  drop_na() %>%
  mutate(phenotype = PA + LS) %>%
  select(hhid, pn, phenotype)

fwrite(SWB, "input/HRS/SWB.pheno")
