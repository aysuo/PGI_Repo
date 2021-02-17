#----------------------------------------------------------------------------------#
# Charts incremental r2 of phenotypes using multitrait scores
# Date: 03/03/2020
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

HRS_path <- "output/prediction/HRS_phenotypes_r2_multi.txt"
WLS_path <- "output/prediction/WLS_phenotypes_r2_multi.txt"

HRS_data <- fread(HRS_path) %>% mutate(Dataset = "HRS")
#
## temporarily drop HRS self-health and asthmaeczrhi
#HRS_data <- filter(
#  HRS_data,
#  phenotype != "SELF_RATED_HEALTH" &
#  phenotype != "AsthmaEczRhi"
#)

WLS_data <- fread(WLS_path) %>% mutate(Dataset = "WLS")

data <- rbind(HRS_data, WLS_data)


########################################################
##################### Wrangle data #####################
########################################################

#unused_phenotypes <- c(
#  "AFB_excl_Barban",
#  "AGE_FIRST_MENSES_excl_Day",
#  "AsthmaEczRhi_excl_Ferreira",
#  "BMI_excl_GIANT", # now have BMI
#  "BMI_test", # now renamed to BMI
#  "EXTRAVERSION_excl_GPC",
#  "HEIGHT_excl_GIANT",
#  "INTELLIGENCE_excl_COGENT",
#  "NEARSIGHTED",
#  "NEUR_excl_deMoor",
#  #"NEBWOM_UKBonly", # now use this for HRS
#  "SWBPA",
#  "SWBLS"
#)
#
## remove NEBwom_UKBonly for WLS only
#data <- data %>%
#  filter(!(phenotype %in% unused_phenotypes))

data <- data %>%
  select(phenotype, r2_inc, r2_inc_lower, r2_inc_upper, Dataset) %>%
  mutate(phenotype = recode(phenotype, "ADHD" = "ADHD_(700k)",
                                       "ADVENTURE" = "Adventurousness_(810k)",
                                       "AFB" = "Age_first_birth_(590k)",
                                       "AGEFIRSTMENSES" = "Age_first_menses_(370k)",
                                       "ALLERGYCAT" = "Cat_allergy_(380k)",
                                       "ALLERGYDUST" = "Dust_allergy_(420k)",
                                       "ALLERGYPOLLEN" = "Pollen_allergy_(330k)",
                                       "ASTECZRHI" = "Asthma,_Eczema and_Rhinitis_(300k)",
                                       "ASTHMA" = "Asthma_(630k)",
                                       #"BMI_excl_GIANT" = "BMI_excluding_GIANT_",
                                       #"BMI_test" = "BMI_test__",
                                       #"BMI" = "BMI_(k)__",
                                       #"CPD" = "Cigarettes_per day_(k)",
                                       "DEP" = "Depression_(900k)",
                                       "EA" = "EA_(1310k)",
                                       #"EVERSMOKE" = "Ever_smoker_(k)",
                                       "EXTRAVERSION" = "Extraversion_(120k)",
                                       "FAMSAT" = "Family_satisfaction_(340k)",
                                       "FINSAT" = "Financial_satisfaction_(500k)",
                                       "HAYFEVER" = "Hay_fever_(600k)",
                                       #"HEIGHT_excl_GIANT" = "Height ex-GIANT",
                                       #"HEIGHT" = "Height_(k)",
                                       "INTELLIGENCE" = "Intelligence_(340k)",
                                       "LEFTOUT" = "Left out_of social_activity_(810k)",
                                       "LONELY" = "Loneliness_(1230k)",
                                       #"MIGRAINE" = "Migraine_(k)",
                                       "NEBmen" = "Number_ever born_(men)_(640k)",
                                       "NEBmen_excl_BarbanNEB" = "Number_ever born_(men excl_Barban NEB)_(450k)",
                                       "NEBwomen" = "Number_ever born_(women)_(450k)",
                                       "NEBwomen_excl_BarbanNEB" = "Number_ever born_(women, excl_Barban NEB)_(360k)",
                                       #"NEUR_excl_deMoor" = "Neuroticism ex-deMoor",
                                       "NEUR" = "Neuroticism_(520k)",
                                       #"OPEN" = "Openness_(k)",
                                       "RELIGIOSITY" = "Religiosity_(880k)",
                                       "RISK" = "Risk_tolerance_(1540k)",
                                       "SELFHEALTH" = "Self-rated_health_(1030k)",
                                       "SWB" = "Subjective_well-being_(1090k)",
                                       "WORKSAT" = "Work_satisfaction_(560k)"),
         phenotype_n = phenotype) %>%
  arrange(phenotype)

########################################################
########### Save data for split-chart script ###########
########################################################

fwrite(data, paste0("output/prediction/HRSandWLS_phenotypes_r2_multi_cleaned.txt"))


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

plot1 <- ggplot(data_plot1, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
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
    "output/prediction/incr2_multi_split_left_",
    Sys.Date(),
    ".png"
  ),
  plot = plot1,
  width = 16,
  height = 6
)

plot2 <- ggplot(data_plot2, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(hjust=0.7),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
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
    "output/prediction/incr2_multi_split_right_",
    Sys.Date(),
    ".png"
  ),
  plot = plot2,
  width = 16,
  height = 6
)


#plot_left <- ggplot(data_left, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
#  theme_minimal() +
#  theme(panel.grid.major.x = element_blank(),
#        axis.title = element_text(size = 17, colour="grey4"),
#        axis.title.x = element_text(hjust=0.65),
#        axis.text = element_text(size = 10.5, colour="grey12"),
#        legend.position = "none") +
#  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
#  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
#  scale_x_discrete(labels = add.line(sort(unique(data_left$phenotype_n)))) +
#  scale_y_continuous(breaks = seq(0,15,2),
#                     limits = c(0, 12.6),
#                     expand = c(0, 0),
#                     labels = function(x) paste0(x, "%")) +
#  xlab("Phenotype") +
#  ylab(expression(paste("Incremental ", italic(R)^2))) +
#  scale_fill_manual(values = c("#1b9e77", "#7570b3"))
#
#plot_right <- ggplot(data_right, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
#  theme_minimal() +
#  theme(panel.grid.major.x = element_blank(),
#        axis.title = element_text(size = 17, colour="grey4"),
#        axis.text = element_text(size = 10.5, colour="grey12")) +
#  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
#  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
#  scale_x_discrete(labels = add.line(sort(unique(data_right$phenotype_n)))) +
#  scale_y_continuous(breaks = seq(0,30,5),
#                     limits = c(0, 31.5),
#                     expand = c(0, 0),
#                     labels = function(x) paste0(x, "%")) +
#  xlab("") +
#  ylab("") +
#  scale_fill_manual(values = c("#1b9e77", "#7570b3"))
#
#grob <- arrangeGrob(
#    plot_left,
#    plot_right,
#
#    nrow=1,
#    widths=c(10, 2.4)
#  )
#
#grid.arrange(grob)
#
#ggsave(
#  file = paste0(
#    "output/prediction/incr2_multi_",
#    Sys.Date(),
#    ".png"
#  ),
#  plot = grob,
#  width = 22,
#  height = 8
#)

plot <- ggplot(data, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(hjust=0.65),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
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
    "output/prediction/incr2_multi_",
    Sys.Date(),
    ".png"
  ),
  plot = plot,
  width = 22,
  height = 8
)
