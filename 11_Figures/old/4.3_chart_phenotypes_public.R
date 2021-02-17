#----------------------------------------------------------------------------------#
# Compares incremental r2 of phenotypes using public and single scores
# Date: 04/22/2020
# Author: Joel Becker

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "dplyr", "ggplot2", "gridExtra", "grid", "gtable")
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
WLS_public_data <- fread(WLS_public_path) %>% mutate(Dataset = "WLS")

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

public_data <- public_data %>%
  rename(
    r2_no_score_public = r2_no_score,
    r2_with_score_public = r2_with_score,
    r2_inc_public = r2_inc,
    r2_inc_lower_public = r2_inc_lower,
    r2_inc_upper_public = r2_inc_upper,
    N_public = N
  )

# changes phenotype names to match
data <- single_data %>%
  mutate(phenotype = recode(phenotype, "NEBwomen" = "NEBWOMEN"))

data <- data %>%
  full_join(public_data, by = c("phenotype", "Dataset")) %>%
  mutate(r2_inc_public = case_when(is.na(r2_inc_public) ~ 0, TRUE ~ r2_inc_public))

# select only desired phenotypes
unused_phenotypes <- c(
  "NEARSIGHTED"
)
data <- data %>%
  filter(!(phenotype %in% unused_phenotypes))

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

########################################################
##################### Wrangle data #####################
########################################################

data <- data %>%
  #select(phenotype, r2_inc, r2_inc_lower, r2_inc_upper, Dataset) %>%
  mutate(
    phenotype = recode(
      phenotype,
      "ACTIVITY" = "Physical_Activity_(90k)",
      "ADHD" = "ADHD_(60k)",
      #"ADVENTURE" = "Adventurousness_(810k)",
      "AFB" = "Age_first_birth_(240k)",
      "ASTECZRHI" = "Asthma,_Eczema and_Rhinitis_(240k)",
      "ASTHMA" = "Asthma_(360k)",
      #"AUDIT" = "AUDIT_(150k)",
      "BMI" = "BMI_(800k)__",
      "CP" = "Intelligence_(270k)",
      "CPD" = "Cigarettes_per day_(260k)",
      "DEP" = "Depression_(500k)",
      "DPW" = "Drinks_Per Week_(540k)",
      "EA" = "EA_(770k)",
      "EVERSMOKE" = "Ever_smoker_(630k)",
      "EXTRA" = "Extraversion_(60k)",
      "FAMSAT" = "Family_satisfaction_(120k)",
      "FRIENDSAT" = "Friend_satisfaction_(120k)",
      "HAYFEVERb" = "Hay_fever_(360k)",
      "HEIGHT" = "Height_(710k)",
      #"LEFTOUT" = "Left out_of social_activity_(810k)",
      "MENARCHE" = "Age_first_menses_(250k)",
      "MIGRAINE" = "Migraine_(360k)",
      "NEBWOMEN" = "Number_ever born_(women)_(240k, 210k)",
      "NEURO" = "Neuroticism_(380k)",
      "OPEN" = "Openness_(20k)",
      "RELIGATT" = "Religiosity_(360k)",
      "RISK" = "Risk_tolerance_(470k)",
      "SELFHEALTH" = "Self-rated_health_(360k)",
      "SWB" = "Subjective_well-being_(200k)",
      "WORKSAT" = "Work_satisfaction_(80k)"
    ),
    phenotype_n = phenotype
  ) %>%
  arrange(phenotype)

########################################################
########### Save data for split-chart script ###########
########################################################

#fwrite(data, paste0("output/prediction/HRSandWLS_phenotypes_r2_multi_cleaned.txt"))


########################################################
################## Annotation function #################
########################################################

add.line <- function(expression){
  gsub('_', '\n', expression)
}

########################################################
################ Separate data by plots ################
########################################################

# separate out left plot vs right_a plot phenotypes
unique_phenotypes <- length(unique(data$phenotype))
n_phenotypes_plot1 <- (unique_phenotypes + 1) %/% 2

phenotypes_plot1 <- unique(data$phenotype)[1:n_phenotypes_plot1]

data_plot1 <- data %>%
  filter(phenotype %in% phenotypes_plot1)

data_plot2 <- data %>%
  filter(!(phenotype %in% phenotypes_plot1))



########################################################
################ Separate data by scale ################
########################################################

#data_left <- data %>%
#  filter(phenotype != "BMI_excluding_GIANT_" & phenotype != "BMI_(760k)__" & phenotype != "EA_(1050k)" & phenotype != "Height_(700k)")
#
#data_right <- data %>%
#  filter(phenotype == "BMI_excluding_GIANT_" | phenotype == "BMI_(760k)__" | phenotype == "EA_(1050k)" | phenotype == "Height_(700k)")


########################################################
###################### Plot chart ######################
########################################################

plot1 <- ggplot(data_plot1, aes(fill=Dataset, x = phenotype, y = r2_inc_public)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower_public, ymax = r2_inc_upper_public), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data_plot1$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,15,2),
                     limits = c(0, 12.6),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2))) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

ggsave(
  file = paste0(
    "output/prediction/incr2_public_split_left_",
    Sys.Date(),
    ".png"
  ),
  plot = plot1,
  width = 16,
  height = 6
)

plot2 <- ggplot(data_plot2, aes(fill=Dataset, x = phenotype, y = r2_inc_public)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(hjust = 0.7),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower_public, ymax = r2_inc_upper_public), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data_plot2$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,15,2),
                     limits = c(0, 12.6),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2))) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

ggsave(
  file = paste0(
    "output/prediction/incr2_public_split_right_",
    Sys.Date(),
    ".png"
  ),
  plot = plot2,
  width = 16,
  height = 6
)

plot <- ggplot(data, aes(fill=Dataset, x = phenotype, y = r2_inc_public)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(hjust=0.65),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower_public, ymax = r2_inc_upper_public), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,15,2),
                     limits = c(0, 12.6),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2))) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

ggsave(
  file = paste0(
    "output/prediction/incr2_public_",
    Sys.Date(),
    ".png"
  ),
  plot = plot,
  width = 22,
  height = 8
)
