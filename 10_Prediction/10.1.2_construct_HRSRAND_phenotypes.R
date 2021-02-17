#!/usr/bin/env Rscript

#----------------------------------------------------------------------------------#
# Constructs phenotypes from HRS and RAND data
# Date: 05/01/2020
# Author: Joel Becker

# Notes:
#  RAND phenotypes constructed first, then HRS phenotypes
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "Rmpfr", "sjmisc")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = T)

# source and set directory
# source("directory_paths.R")
setwd("/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/derived_data/10_Prediction/input/HRS2")
inputDataDir="/disk/genetics/PGS/Aysu/PGS_Repo_pipeline/original_data/prediction_phenotypes/HRS/"

########################################################
################ Residualising functions ###############
########################################################

residualise <- function(data, age_residualise=T, neb=F) {
  # residualise phenotypes within waves and standardise

  # residualise within-wave, record missings
  if (age_residualise==T) { # if residualising on age
    reg <- summary(lm(pheno ~ age + age2 + male + male_age + male_age2, data))
  } else if (neb==F) { # else if residualising on cohort
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

residualise.average.save <- function(data, average=T, age_residualise=T, name, neb=F, rand=F) {
  # call residualising function, optionally average residualise, save output

  # residualise within wave?
  if (average==T) { # yes

    # get data column names for new df
    df <- data[F,]

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
    if (rand == F) {
      df <- df %>%
        group_by(HHID, PN) %>%
        summarise(phenotype = mean(std_resid, na.rm=T)) %>%
        ungroup() %>%
        filter(!is.na(phenotype))
    } else if (rand == T) {
      df <- df %>%
        group_by(hhidpn) %>%
        summarise(phenotype = mean(std_resid, na.rm=T)) %>%
        ungroup() %>%
        filter(!is.na(phenotype))
    }

  } else if (average==F) { # no

    # residualise (ignoring waves)
    data <- residualise(data, age_residualise = age_residualise)

    # select appropriate columns
    if (rand == F) {
      df <- data %>% select(HHID, PN, phenotype=std_resid)
    }

    if (neb==T) { # if NEB, split by sex
      nebmen <- data %>%
        filter(male == 1) %>%
        residualise(., age_residualise = age_residualise, neb=T) %>%
        select(hhidpn, phenotype = std_resid)
      fwrite(nebmen, paste0("NEBmen.pheno"))
      nebwom <- data %>%
        filter(male == 0) %>%
        residualise(., age_residualise = age_residualise, neb=T) %>%
        select(hhidpn, phenotype=std_resid)
      fwrite(nebwom, paste0("NEBwomen.pheno"))
    }

    if (rand == T) {
      df <- data %>% select(hhidpn, phenotype = std_resid)
    }

  }

  # save data
  fwrite(df, paste0(name, ".pheno"))

}


########################################################
#################### Load RAND data ####################
########################################################

RAND_path <- paste0(inputDataDir,"/randhrs1992_2016v2.dta")
RAND_data <- read.dta(RAND_path)


########################################################
############### Construct phenotype: BMI ###############
########################################################

# transform data for residualising
BMI <- RAND_data %>%
  select(hhidpn, contains("agem_m"), contains("gender"), contains("bmi")) %>%
  select(hhidpn, starts_with("r")) %>%
  to_long(keys = "wave", values = c("bmi", "age"),
                                  c(paste0("r", 1:13, "bmi")),
                                  c(paste0("r", 1:13, "agem_m"))) %>%
  mutate(age2 = age^2,
         male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
         male_age = male * age,
         male_age2 = male * age2,
         wave = gsub("r", "", gsub("bmi", "", wave))) %>%
  select(-c(ragender, respagem_m)) %>%
  select(-contains("pmbmi")) %>%
  rename(pheno=bmi) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=BMI, average=T, age_residualise=T, rand = T, name="BMI")

########################################################
############### Construct phenotype: DPW ###############
########################################################

# transform data for residualising
DPW <- RAND_data %>%
  select(
    hhidpn,
    contains("agem_m"),
    contains("gender"),
    contains("drinkd"),
    contains("drinkn")
  ) %>%
  select(hhidpn, starts_with("r"), -"r1agem_m", -"r2agem_m") %>%
  to_long(keys = "wave", values = c("DrinksPerDay", "DaysDrinkingPerWeek", "age"),
                                  c(paste0("r", 3:13, "drinkn")),
                                  c(paste0("r", 3:13, "drinkd")),
                                  c(paste0("r", 3:13, "agem_m"))) %>%
  mutate(
    age2 = age^2,
    male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
    male_age = male * age,
    male_age2 = male * age2,
    wave = gsub("r", "", gsub("dpw", "", wave)),
    pheno = DrinksPerDay * DaysDrinkingPerWeek
  ) %>%
  select(-c(ragender, respagem_m, DrinksPerDay, DaysDrinkingPerWeek)) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=DPW, average=T, age_residualise=T, rand = T, name="DPW")


########################################################
## Construct phenotype: delaydisc (financial horizon) ##
########################################################

# transform data for residualising
delaydisc <- RAND_data %>%
  select(hhidpn,
         contains("agem_m"),
         contains("gender"),
         contains("finpln")) %>%
  select(hhidpn, starts_with("r")) %>%
  select(-c(r2agem_m, r3agem_m, r9agem_m, r10agem_m, respagem_m)) %>%
  to_long(keys = "wave", values = c("finpln", "age"),
                                  c(paste0("r", c(1, 4:8, 11:13), "finpln")),
                                  c(paste0("r", c(1, 4:8, 11:13), "agem_m"))) %>%
  mutate(
    age2 = age^2,
    male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
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

# transform data for residualising
height <- RAND_data %>%
  select(hhidpn,
         contains("agem_m"),
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("height")) %>%
  select(hhidpn, starts_with("r")) %>%
  to_long(keys = "wave", values = c("height", "age"),
                                  c(paste0("r", 1:13, "height")),
                                  c(paste0("r", 1:13, "agem_m"))) %>%
  mutate(height,wave = gsub("r", "", gsub("height", "", wave)),
         male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
         dob = rabyear + (rabmonth / 12),
         dob2 = dob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2) %>%
  select(-c(ragender, respagem_m)) %>%
  group_by(hhidpn) %>%
  mutate(pheno = mean(height, na.rm=T), rn=row_number()) %>%
  ungroup() %>%
  filter(rn==1) %>%
  select(hhidpn, dob, dob2, male, male_dob, male_dob2, pheno) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=height, average=F, age_residualise=F, rand = T, name="HEIGHT")


########################################################
############ Construct phenotype: eversmoke ############
########################################################

# transform data for residualising
smokev <- RAND_data %>%
  select(hhidpn,
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("smokev")) %>%
  select(hhidpn, starts_with("r")) %>%
  mutate(r1smokev = case_when(r1smokev == "1.yes" ~ 1, r1smokev == "0.no" ~ 0),
         r2smokev = case_when(r2smokev == "1.yes" ~ 1, r2smokev == "0.no" ~ 0),
         r3smokev = case_when(r3smokev == "1.yes" ~ 1, r3smokev == "0.no" ~ 0),
         r4smokev = case_when(r4smokev == "1.yes" ~ 1, r4smokev == "0.no" ~ 0),
         r5smokev = case_when(r5smokev == "1.yes" ~ 1, r5smokev == "0.no" ~ 0),
         r6smokev = case_when(r6smokev == "1.yes" ~ 1, r6smokev == "0.no" ~ 0),
         r7smokev = case_when(r7smokev == "1.yes" ~ 1, r7smokev == "0.no" ~ 0),
         r8smokev = case_when(r8smokev == "1.yes" ~ 1, r8smokev == "0.no" ~ 0),
         r9smokev = case_when(r9smokev == "1.yes" ~ 1, r9smokev == "0.no" ~ 0),
         r10smokev = case_when(r10smokev == "1.yes" ~ 1, r10smokev == "0.no" ~ 0),
         r11smokev = case_when(r11smokev == "1.yes" ~ 1, r11smokev == "0.no" ~ 0),
         r12smokev = case_when(r12smokev == "1.yes" ~ 1, r12smokev == "0.no" ~ 0),
         r13smokev = case_when(r13smokev == "1.yes" ~ 1, r13smokev == "0.no" ~ 0)) %>%
  rowwise() %>%
  mutate(smokev = sum(r1smokev, r2smokev, r3smokev, r4smokev,
                      r5smokev, r6smokev, r7smokev, r8smokev,
                      r9smokev, r10smokev, r11smokev, r12smokev,
                      r13smokev, na.rm=T)) %>%
  ungroup() %>%
  mutate(pheno = case_when(smokev >= 1 ~ 1, smokev == 0 ~ 0),
         male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
         dob = rabyear + (rabmonth / 12),
         dob2 = dob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2) %>%
  select(hhidpn, dob, dob2, male, male_dob, male_dob2, pheno) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=smokev, average=F, age_residualise=F, rand = T, name="EVERSMOKE")


########################################################
################ Construct phenotype: EA ###############
########################################################

# transform data for residualising
EA <- RAND_data %>%
  select(hhidpn,
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("raedyrs")) %>%
  select(hhidpn, starts_with("r")) %>%
  mutate(male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
         dob = rabyear + (rabmonth / 12),
         dob2 = dob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2) %>%
  select(hhidpn, dob, dob2, male, male_dob, male_dob2, pheno=raedyrs) %>%
  filter(pheno != "NA") %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=EA, average=F, age_residualise=F, rand = T, name="EA")


########################################################
############ Construct phenotype: loneliness ###########
########################################################

# transform data for residualising
lonely <- RAND_data %>%
  select(hhidpn, contains("agem_m"), contains("gender"), contains("flone")) %>%
  select(hhidpn,
         starts_with("r"),
         r1flone=r1flonex) %>% # note change in first lonely variable, needed for reshape
  to_long(keys = "wave", values = c("flone", "age"),
                                  c(paste0("r", 1:13, "flone")),
                                  c(paste0("r", 1:13, "agem_m"))) %>%
  mutate(wave = gsub("flone", "", wave),
         age2 = age^2,
         male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(flone=="1.all or almost all"~1,
                           flone=="2.most of the time"~(2/3),
                           flone=="3.some of the time"~(1/3),
                           flone=="4.none or almost none"~0,
                           flone == "1.yes" ~ 1,
                           flone == "0.no" ~ 0),
         wave = gsub("r", "", gsub("depres", "", wave))) %>%
  select(-c(ragender, respagem_m)) %>%
  select(-flone) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=lonely, average=T, age_residualise=T, rand = T, name="LONELY")


########################################################
########## Construct phenotype: NEBmen, NEBwom #########
########################################################

# transform data for residualising
neb <- RAND_data %>%
  select(hhidpn,
         contains("rabyear"),
         contains("rabmonth"),
         contains("gender"),
         contains("raevbrn")) %>%
  select(hhidpn, starts_with("r"), -raevbrnf) %>%
  mutate(male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
         dob = rabyear + (rabmonth / 12),
         dob2 = dob^2,
         male_dob = male * dob,
         male_dob2 = male * dob2) %>%
  select(hhidpn, dob, dob2, male, male_dob, male_dob2, pheno=raevbrn) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=neb, average=F, age_residualise=F, name="NEB", rand = T, neb=T)
# note: NEB argument separates data by sex


########################################################
############### Construct phenotype: risk ##############
########################################################

# transform data for residualising
risk <- RAND_data %>%
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
         male = case_when(ragender == "1.male" ~ 1, ragender == "2.female" ~ 0),
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

path_1992 <- paste0(inputDataDir, "1992_merged_respondent_level.txt")
path_1993 <- paste0(inputDataDir, "1993_merged_respondent_level.txt")
path_1994 <- paste0(inputDataDir, "1994_merged_respondent_level.txt")
path_1995 <- paste0(inputDataDir, "1995_merged_respondent_level.txt")
path_1996 <- paste0(inputDataDir, "1996_merged_respondent_level.txt")
path_1998 <- paste0(inputDataDir, "1998_merged_respondent_level.txt")
path_2000 <- paste0(inputDataDir, "2000_merged_respondent_level.txt")
path_2002 <- paste0(inputDataDir, "2002_merged_respondent_level.txt")
path_2004 <- paste0(inputDataDir, "HRS_2004_data.txt")
path_2006 <- paste0(inputDataDir, "HRS_2006_data.txt")
path_2008 <- paste0(inputDataDir, "HRS_2008_data.txt")
path_2010 <- paste0(inputDataDir, "HRS_2010_data.txt")
path_2012 <- paste0(inputDataDir, "HRS_2012_data.txt")
path_2014 <- paste0(inputDataDir, "HRS_2014_data.txt")
path_2016 <- paste0(inputDataDir, "HRS_2016_data.txt")
path_tracker <- paste0(inputDataDir, "TRK2016TR_R_stataversion12.dta")


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

# transform data for residualising
activity <- activity_2004 %>%
 full_join(activity_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 full_join(activity_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
 select(HHID, PN,
        pheno1_2004, pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016,
        pheno2_2004, pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016,
        pheno3_2004, pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016,
        sex,
        mob,
        yob,
        moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
        yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
 gather(key="wave", value="value",
        paste0("pheno1_", seq(2004, 2016, 2)),
        paste0("pheno2_", seq(2004, 2016, 2)),
        paste0("pheno3_", seq(2004, 2016, 2)),
        paste0("moi_",    seq(2004, 2016, 2)),
        paste0("yoi_",    seq(2004, 2016, 2))) %>%
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

# transform data for residualising
asthma <- asthma_2008 %>%
  full_join(asthma_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(asthma_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2008, 2016, 2)),
         paste0("moi_",   seq(2008, 2016, 2)),
         paste0("yoi_",   seq(2008, 2016, 2))) %>%
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

# residualise, average, save
residualise.average.save(data=asthma, average=F, age_residualise=F, name="ASTHMA")


########################################################
############# Construct phenotype: hayfever ############
########################################################

# select and rename relevant variables
hayfever_2008 <- data_2008 %>%
  select(HHID, PN,
         pheno_2008=LB107,
         sex, mob, yob, moi_2008, yoi_2008)

hayfever_2010 <- data_2010 %>%
  select(HHID, PN,
         pheno_2010=MB107,
         sex, mob, yob, moi_2010, yoi_2010)

hayfever_2012 <- data_2012 %>%
  select(HHID, PN,
         pheno_2012=NB107,
         sex, mob, yob, moi_2012, yoi_2012)

hayfever_2014 <- data_2014 %>%
  select(HHID, PN,
         pheno_2014=OB107,
         sex, mob, yob, moi_2014, yoi_2014)

hayfever_2016 <- data_2016 %>%
  select(HHID, PN,
         pheno_2016=PB107,
         sex, mob, yob, moi_2016, yoi_2016)

# transform data for residualising
hayfever <- hayfever_2008 %>%
  full_join(hayfever_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(hayfever_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(hayfever_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(hayfever_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2008, 2016, 2)),
         paste0("moi_",   seq(2008, 2016, 2)),
         paste0("yoi_",   seq(2008, 2016, 2))) %>%
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

# residualise, average, save
#residualise.average.save(data=hayfever, average=F, age_residualise=F, name="HAYFEVER")


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

# transform data for residualising
audit <- audit_2002 %>%
  full_join(audit_2004, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2006, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(audit_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno1_2002, pheno1_2004, pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016,
         pheno2_2002, pheno2_2004, pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016,
         pheno3_2002, pheno3_2004, pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016,
         pheno4_2002, pheno4_2004, pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016,
         pheno5_2002, pheno5_2004, pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016,
         pheno6_2002, pheno6_2004, pheno6_2006, pheno6_2008, pheno6_2010, pheno6_2012, pheno6_2014, pheno6_2016,
         pheno7_2002, pheno7_2004, pheno7_2006, pheno7_2008, pheno7_2010, pheno7_2012, pheno7_2014, pheno7_2016,
         pheno8_2002, pheno8_2004, pheno8_2006, pheno8_2008, pheno8_2010, pheno8_2012, pheno8_2014, pheno8_2016,
         pheno9_2002, pheno9_2004, pheno9_2006, pheno9_2008, pheno9_2010, pheno9_2012, pheno9_2014, pheno9_2016,
         sex,
         mob,
         yob,
         moi_2002, moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2002, yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno1_", seq(2002, 2016, 2)),
         paste0("pheno2_", seq(2002, 2016, 2)),
         paste0("pheno3_", seq(2002, 2016, 2)),
         paste0("pheno4_", seq(2002, 2016, 2)),
         paste0("pheno5_", seq(2002, 2016, 2)),
         paste0("pheno6_", seq(2002, 2016, 2)),
         paste0("pheno7_", seq(2002, 2016, 2)),
         paste0("pheno8_", seq(2002, 2016, 2)),
         paste0("pheno9_", seq(2002, 2016, 2)),
         paste0("moi_",    seq(2002, 2016, 2)),
         paste0("yoi_",    seq(2002, 2016, 2))) %>%
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
  select(HHID, PN,
         pheno_1992, pheno_1993, pheno_1994, pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_1992, moi_1993, moi_1994, moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_1992, yoi_1993, yoi_1994, yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1992:1996, seq(1998, 2016, 2))),
         paste0("moi_",   c(1992:1996, seq(1998, 2016, 2))),
         paste0("yoi_",   c(1992:1996, seq(1998, 2016, 2)))) %>%
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
  select(HHID, PN,

         pheno_1992, pheno_1993, pheno_1994, pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016,

         sex, mob, yob,

         moi_1992, moi_1993, moi_1994, moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016,

         yoi_1992, yoi_1993, yoi_1994, yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1992:1995, seq(1996, 2016, 2))),
         paste0("moi_",   c(1992:1995, seq(1996, 2016, 2))),
         paste0("yoi_",   c(1992:1995, seq(1996, 2016, 2)))) %>%
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
  select(HHID, PN,

         pheno_1995, pheno_1996, pheno_1998,
         pheno_2000, pheno_2002, pheno_2004, pheno_2006, pheno_2008,
         pheno_2010, pheno_2012, pheno_2014, pheno_2016,

         sex,
         mob,
         yob,

         moi_1995, moi_1996, moi_1998,
         moi_2000, moi_2002, moi_2004, moi_2006, moi_2008,
         moi_2010, moi_2012, moi_2014, moi_2016,

         yoi_1995, yoi_1996, yoi_1998,
         yoi_2000, yoi_2002, yoi_2004, yoi_2006, yoi_2008,
         yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", c(1995, seq(1996, 2016, 2))),
         paste0("moi_",   c(1995, seq(1996, 2016, 2))),
         paste0("yoi_",   c(1995, seq(1996, 2016, 2)))) %>%
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

# residualise, average, save
# residualise.average.save(data=migraine, average=F, age_residualise=F, name="MIGRAINE")


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

# transform data for residualising
neuro <- neuro_2006 %>%
  full_join(neuro_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(neuro_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,

         pheno1_2006, pheno1_2008, pheno1_2010,
         pheno1_2012, pheno1_2014, pheno1_2016,

         pheno2_2006, pheno2_2008, pheno2_2010,
         pheno2_2012, pheno2_2014, pheno2_2016,

         pheno3_2006, pheno3_2008, pheno3_2010,
         pheno3_2012, pheno3_2014, pheno3_2016,

         pheno4_2006, pheno4_2008, pheno4_2010,
         pheno4_2012, pheno4_2014, pheno4_2016,

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
         paste0("moi_",    seq(2006, 2016, 2)),
         paste0("yoi_",    seq(2006, 2016, 2))) %>%
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

# transform data for residualising
open <- open_2006 %>%
  full_join(open_2008, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(open_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno1_2006, pheno1_2008, pheno1_2010, pheno1_2012, pheno1_2014, pheno1_2016,
         pheno2_2006, pheno2_2008, pheno2_2010, pheno2_2012, pheno2_2014, pheno2_2016,
         pheno3_2006, pheno3_2008, pheno3_2010, pheno3_2012, pheno3_2014, pheno3_2016,
         pheno4_2006, pheno4_2008, pheno4_2010, pheno4_2012, pheno4_2014, pheno4_2016,
         pheno5_2006, pheno5_2008, pheno5_2010, pheno5_2012, pheno5_2014, pheno5_2016,
         pheno6_2006, pheno6_2008, pheno6_2010, pheno6_2012, pheno6_2014, pheno6_2016,
         pheno7_2006, pheno7_2008, pheno7_2010, pheno7_2012, pheno7_2014, pheno7_2016,
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
         paste0("pheno6_", seq(2006, 2016, 2)),
         paste0("pheno7_", seq(2006, 2016, 2)),
         paste0("moi_",    seq(2006, 2016, 2)),
         paste0("yoi_",    seq(2006, 2016, 2))) %>%
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
  select(HHID, PN,
         pheno_2004, pheno_2006, pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         wave_2004, wave_2006, wave_2008, wave_2010, wave_2012, wave_2014, wave_2016,
         sex,
         mob,
         yob,
         moi_2004, moi_2006, moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2004, yoi_2006, yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
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

# transform data for residualising
selfhealth <- selfhealth_2008 %>%
  full_join(selfhealth_2010, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2012, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2014, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  full_join(selfhealth_2016, by=c("HHID", "PN", "sex", "mob", "yob")) %>%
  select(HHID, PN,
         pheno_2008, pheno_2010, pheno_2012, pheno_2014, pheno_2016,
         sex,
         mob,
         yob,
         moi_2008, moi_2010, moi_2012, moi_2014, moi_2016,
         yoi_2008, yoi_2010, yoi_2012, yoi_2014, yoi_2016) %>%
  gather(key="wave", value="value",
         paste0("pheno_", seq(2008, 2016, 2)),
         paste0("moi_",   seq(2008, 2016, 2)),
         paste0("yoi_",   seq(2008, 2016, 2))) %>%
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

# transform data for residualising
worksat <- worksat_2012 %>%
  mutate(dob = yob + (mob/12),
         doi = yoi_2012 + (moi_2012/12),
         age = doi - dob,
         age2 = age^2,
         male = 2 - sex,
         male_age = male * age,
         male_age2 = male * age2,
         pheno = case_when(pheno_2012==1~0, pheno_2012==2~1, pheno_2012==3~2, pheno_2012==4~3)) %>%
  select(HHID, PN, pheno, age, age2, male, male_age, male_age2) %>%
  drop_na()

# residualise, average, save
residualise.average.save(data=worksat, average=F, age_residualise=T, name="WORKSAT")
