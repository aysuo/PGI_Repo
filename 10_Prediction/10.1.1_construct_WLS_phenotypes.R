#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------#
# Constructs phenotypes from WLS data
# Author: Joel Becker

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "Rmpfr", "sjmisc", "stringr")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# source and set directory
setwd("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction")
inputDataDir="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/WLS/"

########################################################
####################### Load data ######################
########################################################

data <- fread("tmp/WLS_renamed.csv")


########################################################
################ Residualising function ################
########################################################

residualise <- function(data, age_residualise=TRUE, nosex=FALSE) {

  # residualise within-wave, record missings
  if (age_residualise==TRUE) { # if residualising on age
    reg <- summary(lm(pheno ~ age + age2 + male + male_age + male_age2, data))
  } else if (nosex==FALSE) { # else if residualising on cohort
    reg <- summary(lm(pheno ~ dob + dob2 + male + male_dob + male_dob2, data))
  } else {
    reg <- summary(lm(pheno ~ dob + dob2, data))
  }
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

residualise.average.save <- function(data, average=TRUE, age_residualise=TRUE, name, nosex=FALSE, neb=FALSE) {

  # residualise within wave?
  if (average==TRUE) { # yes

    # get data column names for new df
    df <- data[FALSE,]

    # list waves, to run regressions separately
    waves <- sort(unique(data$wave))

    for (i in waves) {
      # residualise within-wave
      data_wave <- filter(data, wave==i)
      data_wave <- residualise(data=data_wave, age_residualise=age_residualise)

      # record this wave
      df <- rbind(df, data_wave)
    }

    # average residuals within-wave
    df <- df %>%
      group_by(id, respondent_type) %>%
      summarise(phenotype = mean(std_resid, na.rm=TRUE)) %>%
      ungroup() %>%
      filter(!is.na(phenotype))

  } else if (average==FALSE) { # no

    # residualise (ignoring waves)
    df <- data %>%
      residualise(., age_residualise = age_residualise, nosex=nosex) %>%
      select(id, respondent_type, phenotype = std_resid)

    if (neb==TRUE) { # if NEB, split by sex
      nebmen <- data %>%
        filter(male == 1) %>%
        residualise(., age_residualise = age_residualise, nosex=neb) %>%
        select(id, respondent_type, phenotype = std_resid)
      nebwom <- data %>%
        filter(male == 0) %>%
        residualise(., age_residualise = age_residualise, nosex=neb) %>%
        select(id, respondent_type, phenotype=std_resid)
      fwrite(nebmen, paste0("input/WLS/NEBmen.pheno"))
      fwrite(nebwom, paste0("input/WLS/NEBwomen.pheno"))
    }

  }

  # save data
  fwrite(df, paste0("input/WLS/", name, ".pheno"))

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

# residualise, average, save
residualise.average.save(data=activity, average=TRUE, age_residualise=TRUE, name="ACTIVITY")


########################################################
############### Construct phenotype: ADHD ##############
########################################################

ADHD <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("ADHD")) %>%
  select(-contains("1957"),
         -contains("1975"),
         -contains("1993"),
         -contains("2004")) %>%
  gather(key = "wave", value = "value",
                               paste0("ADHD_", 2011),
                               paste0("age_", 2011)) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = ADHD) %>% # reverse coded in construction stage
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=ADHD, average=FALSE, age_residualise=TRUE, name="ADHD")


########################################################
######### Construct phenotype: age first birth #########
########################################################

AFB <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("AFB")) %>%
  select(-contains("1957"), -contains("1993"), -contains("2004")) %>%
  mutate(AFB = case_when(!is.na(AFB_2011) ~ AFB_2011, TRUE ~ AFB_1975),
         dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = AFB) %>% # reverse code
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=AFB, average=FALSE, age_residualise=FALSE, name="AFB")


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

# residualise, average, save
residualise.average.save(data=AFM, average=FALSE, age_residualise=FALSE, nosex=TRUE, name="MENARCHE")

########################################################
########## Construct phenotype: agreeableness ##########
########################################################

agree <- data %>%
  select(id_old,
         id,
         respondent_type,
         yob,
         contains("age"),
         male,
         african_american,
         contains("agree")) %>%
  select(-contains("phone"), -contains("nanswered")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key = "wave", value = "value",
                               paste0("agree_", c(1993, 2004, 2011)),
                               paste0("age_",   c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = agree) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=agree, average=TRUE, age_residualise=TRUE, name="AGREE")


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
         paste0("asthma_", c(1993, 2011)), paste0("age_", c(1993, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  group_by(id, respondent_type) %>%
  mutate(asthma = max(asthma, na.rm=TRUE),
         row_n       = row_number()) %>%
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

# residualise, average, save
residualise.average.save(data=asthma, average=FALSE, age_residualise=FALSE, name="ASTHMA")


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

# residualise, average, save
residualise.average.save(data=asthmahayfever, average=FALSE, age_residualise=FALSE, name="ASTECZRHI")



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

# residualise, average, save
residualise.average.save(data=audit, average=T, age_residualise=T, name="AUDIT")


########################################################
############### Construct phenotype: bmi ###############
########################################################

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
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=bmi, average=TRUE, age_residualise=TRUE, name="BMI")


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

# residualise, average, save
residualise.average.save(data=cat, average=FALSE, age_residualise=TRUE, name="ALLERGYCAT")


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

# residualise, average, save
residualise.average.save(data=consc, average=TRUE, age_residualise=TRUE, name="CONSCIENTIOUSNESS")


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
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = copd) %>%
  group_by(id, respondent_type) %>%
  mutate(pheno = max(pheno, na.rm = T)) %>%
  ungroup() %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  filter(pheno %in% 0:1) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=copd, average=TRUE, age_residualise=TRUE, name="COPD")

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

# residualise, average, save
residualise.average.save(data=CPD, average=TRUE, age_residualise=TRUE, name="CPD")


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

# residualise, average, save
residualise.average.save(data=depr, average=TRUE, age_residualise=TRUE, name="DEP")


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

# residualise, average, save
residualise.average.save(data=dpw, average=TRUE, age_residualise=TRUE, name="DPW")


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

# residualise, average, save
residualise.average.save(data=dust, average=FALSE, age_residualise=TRUE, name="ALLERGYDUST")


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

# residualise, average, save
residualise.average.save(data=EA, average=FALSE, age_residualise=FALSE, name="EA")


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

# residualise, average, save
residualise.average.save(data=eversmoke, average=FALSE, age_residualise=FALSE, name="EVERSMOKE")


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

# residualise, average, save
residualise.average.save(data=extra, average=TRUE, age_residualise=TRUE, name="EXTRA")


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

# residualise, average, save
residualise.average.save(data=famsat, average=FALSE, age_residualise=TRUE, name="FAMSAT")


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

# residualise, average, save
residualise.average.save(data=finsat, average=TRUE, age_residualise=TRUE, name="FINSAT")


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

# residualise, average, save
residualise.average.save(data=friendsat1, average=TRUE, age_residualise=TRUE, name="FRIENDSAT1")
residualise.average.save(data=friendsat2, average=TRUE, age_residualise=TRUE, name="FRIENDSAT2")


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

# residualise, average, save
residualise.average.save(data=hayfever, average=FALSE, age_residualise=FALSE, name="HAYFEVER")


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

# residualise, average, save
residualise.average.save(data=height, average=FALSE, age_residualise=FALSE, name="HEIGHT")

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
  drop_na() %>% drop_na()

# residualise, average, save
residualise.average.save(data=intelligence, average=FALSE, age_residualise=FALSE, name="CP")


########################################################
######### Construct phenotype: left out social #########
########################################################

leftoutsocial <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("leftoutsocial"), -contains("phone")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("leftoutsocial_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = leftoutsocial) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=leftoutsocial, average=TRUE, age_residualise=TRUE, name="LEFTOUT")


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

# residualise, average, save
residualise.average.save(data=lonely, average=TRUE, age_residualise=TRUE, name="LONELY")


########################################################
############# Construct phenotype: migraine ############
########################################################

migraine <- data %>%
  select(id_old, id, respondent_type, yob, contains("age"), male, african_american, contains("migraine")) %>%
  select(-contains("1957"), -contains("1975")) %>%
  gather(key="wave", value="value",
         paste0("migraine_", c(1993, 2004, 2011)), paste0("age_", c(1993, 2004, 2011))) %>%
  separate("wave", c("var", "wave")) %>%
  spread("var", "value") %>%
  mutate(age2 = age^2,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = migraine) %>%
  filter(pheno >= 0) %>%
  select(id, respondent_type, wave, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=migraine, average=TRUE, age_residualise=TRUE, name="MIGRAINE")


########################################################
######### Construct phenotype: number ever born ########
########################################################

NEB <- data %>%
  select(id_old, id, respondent_type, yob, male, african_american, contains("NEB")) %>%
  mutate(dob = yob,
         dob2 = yob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2,
         pheno = NEB) %>% # reverse code
  select(id, respondent_type, pheno, dob, dob2, male, male_dob, male_dob2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=NEB, average=FALSE, age_residualise=FALSE, name="NEB", neb=TRUE)


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

# residualise, average, save
residualise.average.save(data=neur, average=TRUE, age_residualise=TRUE, name="NEURO")



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

# residualise, average, save
residualise.average.save(data=open, average=TRUE, age_residualise=TRUE, name="OPEN")


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

# residualise, average, save
residualise.average.save(data=pollen, average=FALSE, age_residualise=TRUE, name="ALLERGYPOLLEN")


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
  drop_na() %>% drop_na()

# residualise, average, save
residualise.average.save(data=relig, average=TRUE, age_residualise=TRUE, name="RELIGATT")


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

# residualise, average, save
residualise.average.save(data=risk, average=TRUE, age_residualise=TRUE, name="RISK")
residualise.average.save(data=risk_loss, average=TRUE, age_residualise=TRUE, name="RISKLOSS")


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

# residualise, average, save
residualise.average.save(data=selfhealth, average=TRUE, age_residualise=TRUE, name="SELFHEALTH")


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

# residualise, average, save
residualise.average.save(data=SWB, average=TRUE, age_residualise=TRUE, name="SWB")


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

# residualise, average, save
residualise.average.save(data=worksat, average=TRUE, age_residualise=TRUE, name="WORKSAT")
