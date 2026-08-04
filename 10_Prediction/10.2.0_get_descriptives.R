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

cohort=args[1]
pheno_file=args[2]
pc_dir=args[3]
crosswalk=args[4]

if ( cohort == "HRS" ){
    birthyearVar <-"rabyear"
    sexVar <- "ragender"
    idVar <- "hhidpn"
} else if ( cohort == "WLS" ){
    birthyearVar <-"yob"
    sexVar <- "male"
    idVar <- "id"
} else if (cohort == "UKB3" ){
    birthyearVar <- "BYEAR"
    sexVar <- "SEX"
    idVar <- "n_eid"
}


if ( grepl(".csv", pheno_file , fixed=TRUE)) {
    pheno_data <- read.csv(pheno_file)
} else if (grepl(".dta", pheno_file , fixed=TRUE) ) {
    pheno_data <- read.dta13(pheno_file)
}

pheno_data <- pheno_data %>%
                rename(ID := !!idVar, 
                        byear := !!birthyearVar,
                        sex := !!sexVar) %>%
                mutate(byear = as.numeric(byear)) %>%
                select(any_of(c("ID","sex","byear","partition"))) 


if ( cohort != "UKB3" ) {
    PCs_data <- fread(paste0(pc_dir,"/",cohort,"_PCs.eigenvec")) %>%
                rename(IID = V2, pc1=V3) %>%
                select("IID","pc1")
}


if ( cohort == "HRS" ) {
    if ( grepl(".csv", crosswalk , fixed=TRUE)) {
        crosswalk_data <- read.csv(crosswalk)
    } else if (grepl(".dta", crosswalk , fixed=TRUE) ) {
        crosswalk_data <- read.dta(crosswalk)
    }
    crosswalk_data <-   crosswalk_data %>%
                        mutate(IID = as.numeric(SUBJID),
                            HHID = as.numeric(HHID),
                            PN = as.numeric(PN)) %>%
                            mutate(ID = (1000 * HHID) + PN)

    pheno_data <- inner_join(pheno_data, crosswalk_data, by="ID")
}

if ( cohort == "WLS" ) {
    pheno_data <- pheno_data %>% 
                    mutate(byear = as.numeric(byear)+1900,
                    IID = paste0(ID,"_",ID))
}

if ( cohort != "UKB3"){
    data <- inner_join(pheno_data, PCs_data, by="IID") %>%
            drop_na(pc1,byear,sex) %>%
            select("IID","byear","sex")
}

if ( cohort == "UKB3" ) {
    data <- pheno_data %>%
            mutate(byear = byear*10+1900) %>%
            filter(partition==3) %>%
            drop_na(byear,sex)
}

describe(data)