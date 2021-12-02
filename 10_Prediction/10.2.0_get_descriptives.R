#!/usr/bin/env Rscript


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "foreign", "dplyr", "tidyr", "Rmpfr", "sjmisc","psych","readstata13")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = T)

args=commandArgs(trailingOnly=TRUE)

########################################################

#### HRS2 #####
cohort="HRS2"
RAND_path <-"original_data/prediction_phenotypes/HRS/randhrs1992_2016v2.dta"
pheno_data <- read.dta(RAND_path)
pheno_data$age <- 2021-pheno_data$rabyear

# PCs
PCs_path <- paste0("derived_data/8_PCs/",cohort,"/",cohort,"_PCs.eigenvec")
    PCs_oldnames <- paste0("V", 3:22)
    PCs_newnames <- paste0("pc", 1:20)
    PCs_data <- fread(PCs_path) %>%
    rename(IID = V2) %>%
    rename_at(vars(PCs_oldnames), ~ PCs_newnames)

# Scores-phenos crosswalk
score_pheno_crosswalk_path <- "original_data/prediction_phenotypes/HRS/HRS_GENOTYPEV2_XREF.dta"
score_pheno_crosswalk_data <- read.dta(score_pheno_crosswalk_path) %>%
    mutate(IID = as.numeric(LOCAL_ID),
    HHID = as.numeric(HHID),
    PN = as.numeric(PN))
score_pheno_crosswalk_data <- score_pheno_crosswalk_data %>% mutate(hhidpn = (1000 * HHID) + PN)

# Merge
data <- inner_join(pheno_data, score_pheno_crosswalk_data, by="hhidpn") %>%
    inner_join(PCs_data, by="IID") %>%
        select(c("IID","age","ragender","pc1")) %>%
        drop_na(pc1,age,ragender)

         
describe(data)

########################################################

#### WLS ####
cohort="WLS"
pheno_data <- fread("derived_data/10_Prediction/tmp/WLS_renamed.csv")
pheno_data$yob <- 1900+pheno_data$yob
pheno_data$age2021 <- 2021-pheno_data$yob

PCs_path <- paste0("derived_data/8_PCs/",cohort,"/",cohort,"_PCs.eigenvec")
    PCs_oldnames <- paste0("V", 3:22)
    PCs_newnames <- paste0("pc", 1:20)
    PCs_data <- fread(PCs_path) %>%
    rename(id = V2) %>%
    rename_at(vars(PCs_oldnames), ~ PCs_newnames)

data <- inner_join(pheno_data, PCs_data, by="id") %>%
        select(c("id","age2021","male","pc1")) %>%
        drop_na(pc1,age2021,male)

describe(data)

########################################################

#### UKB ####
cohort="UKB"
pheno_data <- read.dta13("derived_data/10_Prediction/tmp/pgs_repo.dta")%>% 
        select(c("IID","BYEAR","SEX","partition"))

pheno_data$age<-2021-(pheno_data$BYEAR*10+1900)

describeBy(pheno_data,pheno_data$partition)