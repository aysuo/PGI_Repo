#----------------------------------------------------------------------------------#
# Compares prediction between public and single-trait scores
# Date: 05/06/2020
# Author: Joel Becker

# Notes:
#  TODO: not sure how NEBPOOLED supposed to be working - fixed in 3.5_!
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
setwd(joel_PGS_dir)


########################################################
####################### Load data ######################
########################################################

HRS_single_path <- "output/prediction/HRS_phenotypes_r2_single.txt"
WLS_single_path <- "output/prediction/WLS_phenotypes_r2_single.txt"

HRS_public_path <- "output/prediction/HRS_phenotypes_r2_public.txt"
WLS_public_path <- "output/prediction/WLS_phenotypes_r2_public.txt"

HRS_single_data <- fread(HRS_single_path) %>% mutate(Dataset = "HRS")
WLS_single_data <- fread(WLS_single_path) %>% mutate(Dataset = "WLS")

HRS_public_data <- fread(HRS_public_path) %>% mutate(Dataset = "HRS")

## remove HRS public NEBWOMEN results, which have overlap
#HRS_public_NEBWOMEN <- filter(HRS_public_data, phenotype == "NEBWOMEN")
#HRS_public_NEBWOMEN$r2_no_score <- NA
#HRS_public_NEBWOMEN$r2_with_score <- NA
#HRS_public_NEBWOMEN$r2_inc <- NA
#HRS_public_NEBWOMEN$r2_inc_lower <- NA
#HRS_public_NEBWOMEN$r2_inc_upper <- NA
#HRS_public_data <- HRS_public_data %>%
#  filter(phenotype != "NEBWOMEN") %>%
#  rbind(., HRS_public_NEBWOMEN)
#
## remove HRS public EVERSMOKE results, which have overlap
#HRS_public_EVERSMOKE <- filter(HRS_public_data, phenotype == "EVERSMOKE")
#HRS_public_EVERSMOKE$r2_no_score <- NA
#HRS_public_EVERSMOKE$r2_with_score <- NA
#HRS_public_EVERSMOKE$r2_inc <- NA
#HRS_public_EVERSMOKE$r2_inc_lower <- NA
#HRS_public_EVERSMOKE$r2_inc_upper <- NA
#HRS_public_data <- HRS_public_data %>%
#  filter(phenotype != "EVERSMOKE") %>%
#  rbind(., HRS_public_EVERSMOKE)
#
## remove HRS public EA results, which have overlap
#HRS_public_EA <- filter(HRS_public_data, phenotype == "EA")
#HRS_public_EA$r2_no_score <- NA
#HRS_public_EA$r2_with_score <- NA
#HRS_public_EA$r2_inc <- NA
#HRS_public_EA$r2_inc_lower <- NA
#HRS_public_EA$r2_inc_upper <- NA
#HRS_public_data <- HRS_public_data %>%
#  filter(phenotype != "EA") %>%
#  rbind(., HRS_public_EA)

WLS_public_data <- fread(WLS_public_path) %>% mutate(Dataset = "WLS")

## remove WLS public EA results, which have overlap
#WLS_public_EA <- filter(WLS_public_data, phenotype == "EA")
#WLS_public_EA$r2_no_score <- NA
#WLS_public_EA$r2_with_score <- NA
#WLS_public_EA$r2_inc <- NA
#WLS_public_EA$r2_inc_lower <- NA
#WLS_public_EA$r2_inc_upper <- NA
#WLS_public_data <- WLS_public_data %>%
#filter(phenotype != "EA") %>%
#rbind(., WLS_public_EA)

single_data <- rbind(HRS_single_data, WLS_single_data)
public_data <- rbind(HRS_public_data, WLS_public_data)

# rename columns for single/public respectively
single_data <- single_data %>%
  rename(
    r2_no_score_single = r2_no_score,
    r2_with_score_single = r2_with_score,
    r2_inc_single = r2_inc,
    r2_inc_lower_single = r2_inc_lower,
    r2_inc_upper_single = r2_inc_upper,
    N_single = N
  )

# append NEBWOM single-scores data to NEBPOOLED public-scores data
#WLS_NEBPOOLED <- single_data %>%
#  filter(Dataset == "WLS" & phenotype == "NEBwomen") %>%
#  mutate(phenotype = "NEBPOOLED")
#
#single_data <- rbind(single_data, WLS_NEBPOOLED)

public_data <- public_data %>%
  rename(
    r2_no_score_public = r2_no_score,
    r2_with_score_public = r2_with_score,
    r2_inc_public = r2_inc,
    r2_inc_lower_public = r2_inc_lower,
    r2_inc_upper_public = r2_inc_upper,
    N_public = N
  ) %>%
  mutate(phenotype = recode(phenotype, "NEBWOMEN" = "NEBwomen"))

# changes phenotype names to match
data <- single_data #%>%
#  mutate(phenotype = recode(phenotype,
#                            #"SELFHEALTH" = "SELF_RATED_HEALTH",
#                            "RELIGATT" = "RELIGIOSITY",
#                            #"NEURO" = "NEUR",
#                            "NEBwomen" = "NEBWOMEN",
#                            "MENARCHE" = "AGEFIRSTMENSES",
#                            "EXTRA" = "EXTRAVERSION",
#                            "DEP" = "DEPRESSION",
#                            "AsthmaEczRhi" = "ASTECZRHI"))

data <- data %>%
  full_join(public_data, by = c("phenotype", "Dataset")) %>%
  mutate(r2_inc_public = case_when(is.na(r2_inc_public) ~ 0, TRUE ~ r2_inc_public))

# select only desired phenotypes
unused_phenotypes <- c(
  "NEARSIGHTED"
)

# remove NEBwom_UKBonly for WLS only
data <- data %>%
  # select only desired phenotypes
  filter(!(phenotype %in% unused_phenotypes)) %>%
  # remove phenotypes where cohorts overlap
  #filter(!(phenotype == "EA"        & Dataset == "HRS")) %>%
  #filter(!(phenotype == "EA"        & Dataset == "WLS")) %>%
  #filter(!(phenotype == "EVERSMOKE" & Dataset == "HRS")) %>%
  #filter(!(phenotype == "SWB"            & Dataset == "HRS")) %>%
  # remove redundant phenotypes
  #filter(!(phenotype == "NEBwomen"  & Dataset == "HRS")) %>%
  # format digits
  mutate(
    r2_inc_single = case_when(
      r2_inc_single != 0 ~ formatC(r2_inc_single, digits = 2, format = "f"),
      r2_inc_single == 0 ~ formatC(r2_inc_single, digits = 0, format = "f")
    ),
    r2_inc_lower_single = case_when(
      r2_inc_lower_single != 0 ~ formatC(r2_inc_lower_single, digits = 2, format = "f"),
      r2_inc_lower_single == 0 ~ formatC(r2_inc_lower_single, digits = 0, format = "f")
    ),
    r2_inc_upper_single = case_when(
      r2_inc_upper_single != 0 ~ formatC(r2_inc_upper_single, digits = 2, format = "f"),
      r2_inc_upper_single == 0 ~ formatC(r2_inc_upper_single, digits = 0, format = "f")
    ),
    r2_inc_public = case_when(
      r2_inc_public != 0 ~ formatC(r2_inc_public, digits = 2, format = "f"),
      r2_inc_public == 0 ~ formatC(r2_inc_public, digits = 0, format = "f")
    ),
    r2_inc_lower_public = case_when(
      r2_inc_lower_public != 0 ~ formatC(r2_inc_lower_public, digits = 2, format = "f"),
      r2_inc_lower_public == 0 ~ formatC(r2_inc_lower_public, digits = 0, format = "f")
    ),
    r2_inc_upper_public = case_when(
      r2_inc_upper_public != 0 ~ formatC(r2_inc_upper_public, digits = 2, format = "f"),
      r2_inc_upper_public == 0 ~ formatC(r2_inc_upper_public, digits = 0, format = "f")
    )
  ) %>%
  arrange(phenotype, Dataset)

select(
  data,
  phenotype,
  Dataset,
  single  = r2_inc_single,
  public  = r2_inc_public
)


########################################################
###################### Save data #######################
########################################################

fwrite(data,
  paste0(getwd(),
    "/output/prediction/public_single_comparison_r2_",
    Sys.Date(),
    ".txt")
  )

# overwrite main file
fwrite(data,
  paste0(getwd(),
    "/output/prediction/public_single_comparison_r2.txt")
  )
