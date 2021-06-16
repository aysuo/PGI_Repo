#----------------------------------------------------------------------------------#
# Generates intermediate data to use for plotting
# Date: 04/16/2019
# Author: Joel Becker, Hariharan Jayashankar

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "readxl", "dplyr", "ggplot2", "gridExtra", "grid", "gtable")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

args=commandArgs(trailingOnly=TRUE)
mainDir=args[1]

inpath <- paste0(mainDir,"/code/11_Figures/")
outpath <- paste0(mainDir,"/derived_data/11_Figures/Prediction/input/") 
excelpath <- paste0(inpath, "Supplementary Tables.xlsx")
sheet <- "3. Prediction results"

data_phenotype <- read_excel(excelpath,
                             sheet = sheet,
                             range = "B4:B65",
                             col_names = "phenotype")

# HRS ---------------------------------------------------------------------
data_HRS <- read_excel(excelpath,
                   sheet = sheet,
                   range = "D4:I65",
                   col_names = c("r2_inc_single",
                                 "r2_inc_multi",
                                 "r2_inc_public",
                                 "r2_diff_single_public",
                                 "r2_diff_multi_public",
                                 "r2_diff_multi_single"))

data_HRS <- cbind(data_phenotype, data_HRS)

data_HRS <- data_HRS %>%
    mutate(r2_inc_single_ci = str_extract(r2_inc_single, "\\([^()]+\\)"),
           r2_inc_multi_ci = str_extract(r2_inc_multi, "\\([^()]+\\)"),
           r2_inc_public_ci = str_extract(r2_inc_public, "\\([^()]+\\)"),
           r2_diff_single_public_ci = str_extract(r2_diff_single_public, "\\([^()]+\\)"),
           r2_diff_multi_public_ci = str_extract(r2_diff_multi_public, "\\([^()]+\\)"),
           r2_diff_multi_single_ci = str_extract(r2_diff_multi_single, "\\([^()]+\\)"),
           
           # lower ci
           r2_inc_single_lower = str_extract(r2_inc_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_lower = str_extract(r2_inc_multi_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_lower = str_extract(r2_inc_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_lower = str_extract(r2_diff_single_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_lower = str_extract(r2_diff_multi_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_lower = str_extract(r2_diff_multi_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           
           # Upper
           r2_inc_single_upper = str_extract(r2_inc_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_upper = str_extract(r2_inc_multi_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_upper = str_extract(r2_inc_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_upper = str_extract(r2_diff_single_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_upper = str_extract(r2_diff_multi_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_upper = str_extract(r2_diff_multi_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric())

# redefine var
data_HRS <-data_HRS %>%
    mutate(r2_inc_single = str_extract(r2_inc_single, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_multi = str_extract(r2_inc_multi, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_public = str_extract(r2_inc_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_single_public = str_extract(r2_diff_single_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_public = str_extract(r2_diff_multi_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_single = str_extract(r2_diff_multi_single, "[-]?\\d+.\\d+") %>% as.numeric())
    
data_HRS <- data_HRS %>% 
    mutate(Dataset = "HRS")


# WLS ---------------------------------------------------------------------
data_WLS <- read_excel(excelpath,
                       sheet = sheet,
                       range = "K4:P65",
                       col_names = c("r2_inc_single",
                                     "r2_inc_multi",
                                     "r2_inc_public",
                                     "r2_diff_single_public",
                                     "r2_diff_multi_public",
                                     "r2_diff_multi_single"))

data_WLS <- cbind(data_phenotype, data_WLS)

data_WLS <- data_WLS %>%
    mutate(r2_inc_single_ci = str_extract(r2_inc_single, "\\([^()]+\\)"),
           r2_inc_multi_ci = str_extract(r2_inc_multi, "\\([^()]+\\)"),
           r2_inc_public_ci = str_extract(r2_inc_public, "\\([^()]+\\)"),
           r2_diff_single_public_ci = str_extract(r2_diff_single_public, "\\([^()]+\\)"),
           r2_diff_multi_public_ci = str_extract(r2_diff_multi_public, "\\([^()]+\\)"),
           r2_diff_multi_single_ci = str_extract(r2_diff_multi_single, "\\([^()]+\\)"),
           
           # lower ci
           r2_inc_single_lower = str_extract(r2_inc_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_lower = str_extract(r2_inc_multi_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_lower = str_extract(r2_inc_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_lower = str_extract(r2_diff_single_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_lower = str_extract(r2_diff_multi_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_lower = str_extract(r2_diff_multi_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           
           # Upper
           r2_inc_single_upper = str_extract(r2_inc_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_upper = str_extract(r2_inc_multi_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_upper = str_extract(r2_inc_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_upper = str_extract(r2_diff_single_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_upper = str_extract(r2_diff_multi_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_upper = str_extract(r2_diff_multi_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric())

# redefine var
data_WLS <-data_WLS %>%
    mutate(r2_inc_single = str_extract(r2_inc_single, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_multi = str_extract(r2_inc_multi, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_public = str_extract(r2_inc_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_single_public = str_extract(r2_diff_single_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_public = str_extract(r2_diff_multi_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_single = str_extract(r2_diff_multi_single, "[-]?\\d+.\\d+") %>% as.numeric())


data_WLS <- data_WLS %>% 
    mutate(Dataset = "WLS")


# Dunedin -----------------------------------------------------------------
data_dunedin <- read_excel(excelpath,
                       sheet = sheet,
                       range = "R4:W65",
                       col_names = c("r2_inc_single",
                                     "r2_inc_multi",
                                     "r2_inc_public",
                                     "r2_diff_single_public",
                                     "r2_diff_multi_public",
                                     "r2_diff_multi_single"))

data_dunedin <- cbind(data_phenotype, data_dunedin)

data_dunedin <- data_dunedin %>%
    mutate(r2_inc_single_ci = str_extract(r2_inc_single, "\\([^()]+\\)"),
           r2_inc_multi_ci = str_extract(r2_inc_multi, "\\([^()]+\\)"),
           r2_inc_public_ci = str_extract(r2_inc_public, "\\([^()]+\\)"),
           r2_diff_single_public_ci = str_extract(r2_diff_single_public, "\\([^()]+\\)"),
           r2_diff_multi_public_ci = str_extract(r2_diff_multi_public, "\\([^()]+\\)"),
           r2_diff_multi_single_ci = str_extract(r2_diff_multi_single, "\\([^()]+\\)"),
           
           # lower ci
           r2_inc_single_lower = str_extract(r2_inc_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_lower = str_extract(r2_inc_multi_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_lower = str_extract(r2_inc_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_lower = str_extract(r2_diff_single_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_lower = str_extract(r2_diff_multi_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_lower = str_extract(r2_diff_multi_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           
           # Upper
           r2_inc_single_upper = str_extract(r2_inc_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_upper = str_extract(r2_inc_multi_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_upper = str_extract(r2_inc_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_upper = str_extract(r2_diff_single_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_upper = str_extract(r2_diff_multi_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_upper = str_extract(r2_diff_multi_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric())

# redefine var
data_dunedin <- data_dunedin %>%
    mutate(r2_inc_single = str_extract(r2_inc_single, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_multi = str_extract(r2_inc_multi, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_public = str_extract(r2_inc_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_single_public = str_extract(r2_diff_single_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_public = str_extract(r2_diff_multi_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_single = str_extract(r2_diff_multi_single, "[-]?\\d+.\\d+") %>% as.numeric())


data_dunedin <- data_dunedin %>% 
    mutate(Dataset = "Dunedin")


# E-risk -----------------------------------------------------------------
data_erisk <- read_excel(excelpath,
                           sheet = sheet,
                           range = "Y4:AD65",
                           col_names = c("r2_inc_single",
                                         "r2_inc_multi",
                                         "r2_inc_public",
                                         "r2_diff_single_public",
                                         "r2_diff_multi_public",
                                         "r2_diff_multi_single"))

data_erisk <- cbind(data_phenotype, data_erisk)

data_erisk <- data_erisk %>%
    mutate(r2_inc_single_ci = str_extract(r2_inc_single, "\\([^()]+\\)"),
           r2_inc_multi_ci = str_extract(r2_inc_multi, "\\([^()]+\\)"),
           r2_inc_public_ci = str_extract(r2_inc_public, "\\([^()]+\\)"),
           r2_diff_single_public_ci = str_extract(r2_diff_single_public, "\\([^()]+\\)"),
           r2_diff_multi_public_ci = str_extract(r2_diff_multi_public, "\\([^()]+\\)"),
           r2_diff_multi_single_ci = str_extract(r2_diff_multi_single, "\\([^()]+\\)"),
           
           # lower ci
           r2_inc_single_lower = str_extract(r2_inc_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_lower = str_extract(r2_inc_multi_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_lower = str_extract(r2_inc_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_lower = str_extract(r2_diff_single_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_lower = str_extract(r2_diff_multi_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_lower = str_extract(r2_diff_multi_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           
           # Upper
           r2_inc_single_upper = str_extract(r2_inc_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_upper = str_extract(r2_inc_multi_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_upper = str_extract(r2_inc_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_upper = str_extract(r2_diff_single_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_upper = str_extract(r2_diff_multi_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_upper = str_extract(r2_diff_multi_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric())

# redefine var
data_erisk <- data_erisk %>%
    mutate(r2_inc_single = str_extract(r2_inc_single, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_multi = str_extract(r2_inc_multi, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_public = str_extract(r2_inc_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_single_public = str_extract(r2_diff_single_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_public = str_extract(r2_diff_multi_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_single = str_extract(r2_diff_multi_single, "[-]?\\d+.\\d+") %>% as.numeric())


data_erisk <- data_erisk %>% 
    mutate(Dataset = "E-Risk")


# UKB ---------------------------------------------------------------------
data_UKB <- read_excel(excelpath,
                         sheet = sheet,
                         range = "AF4:AK65",
                         col_names = c("r2_inc_single",
                                       "r2_inc_multi",
                                       "r2_inc_public",
                                       "r2_diff_single_public",
                                       "r2_diff_multi_public",
                                       "r2_diff_multi_single"))

data_UKB <- cbind(data_phenotype, data_UKB)

data_UKB <- data_UKB %>%
    mutate(r2_inc_single_ci = str_extract(r2_inc_single, "\\([^()]+\\)"),
           r2_inc_multi_ci = str_extract(r2_inc_multi, "\\([^()]+\\)"),
           r2_inc_public_ci = str_extract(r2_inc_public, "\\([^()]+\\)"),
           r2_diff_single_public_ci = str_extract(r2_diff_single_public, "\\([^()]+\\)"),
           r2_diff_multi_public_ci = str_extract(r2_diff_multi_public, "\\([^()]+\\)"),
           r2_diff_multi_single_ci = str_extract(r2_diff_multi_single, "\\([^()]+\\)"),
           
           # lower ci
           r2_inc_single_lower = str_extract(r2_inc_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_lower = str_extract(r2_inc_multi_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_lower = str_extract(r2_inc_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_lower = str_extract(r2_diff_single_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_lower = str_extract(r2_diff_multi_public_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_lower = str_extract(r2_diff_multi_single_ci, "\\([-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           
           # Upper
           r2_inc_single_upper = str_extract(r2_inc_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_multi_upper = str_extract(r2_inc_multi_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_inc_public_upper = str_extract(r2_inc_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_single_public_upper = str_extract(r2_diff_single_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_public_upper = str_extract(r2_diff_multi_public_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric(),
           r2_diff_multi_single_upper = str_extract(r2_diff_multi_single_ci, " [-]?\\d+.\\d+") %>%
               str_sub(2, 6) %>% as.numeric())

# redefine var
data_UKB <- data_UKB %>%
    mutate(r2_inc_single = str_extract(r2_inc_single, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_multi = str_extract(r2_inc_multi, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_inc_public = str_extract(r2_inc_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_single_public = str_extract(r2_diff_single_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_public = str_extract(r2_diff_multi_public, "[-]?\\d+.\\d+") %>% as.numeric(),
           r2_diff_multi_single = str_extract(r2_diff_multi_single, "[-]?\\d+.\\d+") %>% as.numeric())



data_UKB <- data_UKB %>% 
    mutate(Dataset = "UKB3")



# Joining Datasets and Cleaning -------------------------------------------

dfout <- rbind(data_HRS, data_WLS, data_dunedin, data_erisk, data_UKB)

dfout <- dfout %>%
    drop_na(phenotype)

# phenotype labels
# manually put in values here from sample_size_single and sample_size_multi
phenotype_labels <- tribble(
    ~phenotype,    ~phenotype_label,          ~sample_single,       ~sample_multi,
    "Physical Activity" ,   "Physical_Activity",       "357k, 357k,_357k, 357k",               "",
    "Attention Deficit Hyperactivity Disorder (ADHD)" ,       "ADHD",                    "117k, 117k",               "760k, 760k,_760k, 760k, 607k",
    "Adventurousness" ,  "Adventur-_ousness",       "557k, 557k",               "801k, 801k,_801k, 801k, 760k",
    "Age First Birth" ,        "Age First_Birth",         "407k, 407k, 142k",               "368k, 373k,_373k, 373k, 345k",
    "Allergy - Pollen", "Allergy:_Pollen",       "",                   "249k, 249k,_249k, 249k, 214k",
    "Allergy - Dust",   "Allergy:_Dust",          "",                  "355k, 355k,_355k, 355k, 265k",
    "Allergy - Cat",    "Allergy:_Cat",            "",                  "369k, 369k,_369k, 369k, 252k",
    "Asthma" ,     "Asthma",                  "445k, 445k,_445k, 445k, 297k",               "560k, 560k,_560k, 560k, 380k",
    "Asthma/Eczema/Rhinitis" ,  "Asthma/_Eczema/_Rhinitis", "685k, 685k, 433k",              "925k, 925k,_925k, 925k, 674k",
    "Alcohol Misuse" ,      "Alcohol Misuse",          "151k, 151k, 110k",      "303k, 371k,_371k, 371k, 202k",
    "Body Mass Index (BMI)" ,        "BMI",                     "760k, 760k,_760k, 760k, 615k",               "",
    "Cannabis Use",         "Cannabis Use",            "158k",                "",
    "Cognitive Empathy",          "Cognitive Empathy",       "",                                   "56k, 56k,_56k, 56k, 56k",
    "Cognitive Performance" ,         "Cognitive_Performance",   "260k, 260k,_260k, 188k",           "343k, 343k,_343k, 343k, 277k",
    "Delay Disounting",   "Delay Discounting",        "",                   "445k, 445k,_445k, 445k, 413k",
    "Cigarettes per Day" ,        "Cigarettes_per Day",      "254k, 340k,_340k, 207k",           "",
    "COPD",         "COPD",                   "",                "1869k, 1869k,_1869k, 1869k, 1212k",
    "Depressive Symptoms" ,        "Depressive_Symptoms",     "943k, 943k,_943k, 943k, 808k",           "1306k, 1306k,_1306k, 1306k, 1102k",
    "Drinks per Week" ,        "Drinks_per Week",         "385k, 537k,_537k, 537k, 257k",           "460k, 617k,_617k, 617k, 306k",
    "Educational Attainment" ,         "EA",                      "1047k, 1047k,_1047k, 1047k, 984k",          "1295k, 1295k,_1295k, 1295k, 1021k",
    "Ever Smoker" ,  "Ever_Smoker",             "1143k, 1256k,_1256k, 1256k, 994k",           "",
    "Extraversion" ,      "Extraversion",            "122k, 122k,_122k, 122k",           "111k",
    "Life Satisfaction: Family" ,     "Life Satisfaction:_Family", "168k, 168k, 115k",         "290k, 290k,_290k, 290k, 255k",
    "Life Satisfaction: Finance",       "Life Satisfaction:_Finance", "",           "491k, 491k,_491k, 491k, 490k",
    "Life Satisfaction: Friends" ,  "Life Satisfaction:_Friend", "168k, 114k",         "402k, 402k,_402k, 402k, 336k",
    "Hayfever" ,   "Hayfever",                "461k, 461k, 312k",           "284k, 284k,_284k, 284k, 202k",
    "Height" ,     "Height",                  "698k, 698k,_698k, 698k, 550k",           "",
    "Highest Math",        "Highest Math",            "430k", "801k,_801k, 801k, 801k, 775k",
    "Left Out of Social Activity" ,    "Left Out_of Social_Activity", "507k",       "802k, 802k,_802k, 802k, 780k",
    "Loneliness",      "Loneliness",               "",              "1170k, 1170k,_1170k, 1170k, 1138k",
    "Age First Menses (women)" ,   "Age First_Menses",        "329k, 367k",           "320k, 320k,_320k, 320k, 363k",
    "Migraine" ,   "Migraine",                "693k, 558k",           "",
    "Morning Person",     "Morning Person",           "493k, 360k",             "", ###
    "Narcissism",        "Narcissism",              "452k", "",
    "Nearsightedness",        "Nearsightedness",         "367k, 367k,_367k, 367k, 316k",               "", ###
    "Number Ever Born (men)",       "Number_Ever Born_(Men)", "",               "488k, 594k,_594k, 594k, 332k",
    "Number Ever Born (women)" ,   "Number_Ever Born_(women)", "241k, 399k,_399k, 162k",    "464k, 497k,_497k, 497k, 267k",
    "Neuroticism" ,      "Neuroticism",             "485k, 485k,_485k, 485k, 364k",           "480k, 480k,_480k, 480k, 386k",
    "Openness" ,       "Openness",                "76k, 76k, 76k, 76k",            "",
    "Religious Attendance" ,   "Religious_Attendance",    "445k, 445k, 297k",           "793k, 793k,_793k, 793k, 632k",
    "Childhood Reading",     "Childhood_Reading",     "173k, 173k",                "",
    "Risk" ,       "Risk_Tolerance",                   "1427k, 1427k,_1427k, 1257k",           "1753k, 1753k,_1753k, 1753k, 1605k",
    "Self-Rated Health" , "Self-Rated_Health",      "1203k, 1203k,_1203k, 1055k",  "1250k, 1250k,_1250k, 1250k, 1075k",
    "Self-Rated Math Ability",     "Self-Rated_Math Ability", "565k", "663k, 663k,_663k, 663k, 651k",
    "Subjective Well-Being" ,        "Subjective_Well-Being",  "1023k, 1023k,_1023k, 1023k, 969k",           "1619k, 1619k,_1619k, 1619k, 1508k",
    "Life Satisfaction: Work",     "Life Satisfaction:_Work",       "",         "702k",
    "Age Voice Deepened (men)",       "Age Voice_Deepened",   "",      "295k, 295k,_295k, 295k, 301k"
    
)

phenotype_codes <- tibble(
    phenotype_code = c(
        "ACTIVITY" ,
        "ADHD" ,
        "ADVENTURE" ,
        "AFB" ,
        "ALLERGYPOLLEN",
        "ALLERGYDUST",
        "ALLERGYCAT",
        "ASTHMA" ,
        "ASTECZRHI" ,
        "HHAC" ,
        "BMI" ,
        "CANNABIS",
        "CE",
        "CP" ,
        "DELAYDISC",
        "CPD" ,
        "COPD",
        "DEP" ,
        "DPW" ,
        "EA" ,
        "EVERSMOKE" ,
        "EXTRA" ,
        "FAMSAT" ,
        "FINSAT",
        "FRIENDSAT" ,
        "HAYFEVER" ,
        "HEIGHT" ,
        "MATH",
        "LEFTOUT" ,
        "LONELY",
        "MENARCHE" ,
        "MIGRAINE" ,
        "MORNING",
        "NARC",
        "NEAR",
        "NEBmen",
        "NEBwomen" ,
        "NEURO" ,
        "OPEN" ,
        "RELIGATT" ,
        "READING",
        "RISK" ,
        "SELFHEALTH" ,
        "SELFMATH",
        "SWB" ,
        "WORKSAT",
        "VOICE"
    )
)

phenotype_labels <- cbind(phenotype_labels, phenotype_codes)

phenotype_labels <- phenotype_labels %>% 
    rename(phenotype_label2 = phenotype,
           phenotype = phenotype_code)

phenotype_labels %>% fwrite("outdata/phenotype_labels.txt")



dfout <- dfout %>%
    left_join(phenotype_labels %>% rename(phenotype_code = phenotype,
                                          phenotype = phenotype_label2), by = "phenotype")

# making groups
dfout <- dfout %>%
    mutate(pheno_group = case_when(phenotype_code %in% c("BMI", 
                                                      "HEIGHT") ~ "Anthropometric",
                                   phenotype_code %in% c("CP",
                                                         "EA",
                                                      "READING") ~ "Cognition",
                                   phenotype_code %in% c("AFB",
                                                      "MENARCHE",
                                                      "NEBmen",
                                                      "NEBwomen") ~ 
                                       "Fertility and Sexual Development",
                                   phenotype_code %in% c("ADHD",
                                                      "ALLERGYPOLLEN",
                                                      "ALLERGYDUST",
                                                      "ALLERGYCAT",
                                                      "ASTHMA",
                                                      "ASTECZRHI",
                                                      "COPD",
                                                      "CPD",
                                                      "DEP",
                                                      "EVERSMOKE",
                                                      "HAYFEVER",
                                                      "MIGRAINE",
                                                      "SELFHEALTH",
                                                      "DPW",
                                                      "ACTIVITY",
                                                      "HHAC",
                                                      "CANNABIS") ~
                                       "Health and Health Behaviors",
                                   phenotype_code %in% c("EXTRA",
                                                      "LEFTOUT",
                                                      "FAMSAT",
                                                      "FINSAT",
                                                      "FRIENDSAT",
                                                      "WORKSAT",
                                                      "LONELY",
                                                      "NEURO",
                                                      "RELIGATT",
                                                      "RISK",
                                                      "SWB",
                                                      "OPEN",
                                                      "ADVENTURE",
                                                      "MORNING") ~
                                       "Personality and Well-Being"))


cat("Checking if all phenotypes have a phenotype group")
stopifnot(data %>% filter(is.na(pheno_group)) %>% summarize(n()) %>% pull(`n()`) == 0)


# Ordering of data
dfout <- dfout %>%
    mutate(phenotype_n = factor(phenotype_code,
                                levels = c("READING",
                                           "BMI", 
                                           "HEIGHT",
                                           
                                           "CP",
                                           "EA",
                                           
                                           "AFB",
                                           "MENARCHE",
                                           "NEBmen",
                                           "NEBwomen",
                                           
                                           "ADHD",
                                           "HHAC",
                                           "ALLERGYPOLLEN",
                                           "ALLERGYDUST",
                                           "ALLERGYCAT",
                                           "ASTHMA",
                                           "ASTECZRHI",
                                           "CANNABIS",
                                           "COPD",
                                           "CPD",
                                           "DEP",
                                           "DPW",
                                           "EVERSMOKE",
                                           "HAYFEVER",
                                           "MIGRAINE",
                                           "ACTIVITY",
                                           "SELFHEALTH",
                                           
                                           "ADVENTURE",
                                           "EXTRA",
                                           "LEFTOUT",
                                           "FAMSAT",
                                           "FINSAT",
                                           "FRIENDSAT",
                                           "WORKSAT",
                                           "LONELY",
                                           "MORNING",
                                           "NEURO",
                                           "OPEN",
                                           "RELIGATT",
                                           "RISK",
                                           "SWB"
                                )))


dfout <- dfout %>%
    mutate(Dataset = factor(Dataset,
                            levels = c("HRS", "WLS", "Dunedin", "E-Risk", "UKB3")))


dfout <- dfout %>%
    rename(phenotype_label2 = phenotype)

dfout <- dfout %>%
    mutate(phenotype = phenotype_n)

########################################################
########### Save data for split-chart script ###########
########################################################

fwrite(dfout, paste0(outpath, "phenotypes_r2_cleaned.txt"))
saveRDS(dfout, paste0(outpath, "phenotypes_r2_cleaned.R"))
