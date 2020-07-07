#----------------------------------------------------------------------------------#
# Charts incremental r2 of phenotypes
# Date: 04/16/2019
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

HRS_path <- "output/prediction/HRS_phenotypes_r2_single.txt"
WLS_path <- "output/prediction/WLS_phenotypes_r2_single.txt"

HRS_data <- fread(HRS_path) %>% mutate(Dataset = "HRS")
WLS_data <- fread(WLS_path) %>% mutate(Dataset = "WLS")

# temp HRS change column selection
HRS_data <- HRS_data[, colnames(WLS_data)]

data <- rbind(HRS_data, WLS_data)


########################################################
##################### Wrangle data #####################
########################################################

unused_phenotypes <- c(
  "NEARSIGHTED"
)

data <- data %>%
  select(phenotype, r2_inc, r2_inc_lower, r2_inc_upper, Dataset) %>%
  mutate(
    #phenotype = case_when(
    #  phenotype == "NEBwomen" & Dataset == "HRS"
    #),
    phenotype = recode(phenotype,
      "ACTIVITY" = "Physical_Activity_(360k)",
      "ADHD" = "ADHD_(110k)",
      "ADVENTURE" = "Adventurousness_(560k)",
      "AFB" = "Age_First_Birth_(390k)",
      "ASTHMA" = "Asthma_(450k)",
      "ASTECZRHI" = "Asthma,_Eczema and_Rhinitis_(320k)",
      "AUDIT" = "AUDIT_(150k)",
      "BMI" = "BMI_(760k)__",
      "CP" = "Cognitive_Performance_(260k)",
      "CPD" = "Cigarettes_per day_(140k)",
      "DEP" = "Depression_(740k)",
      "DPW" = "Drinks_Per Week_(540k)",
      "EA" = "EA_(1050k)",
      "EVERSMOKE" = "Ever_Smoker_(520k)",
      "EXTRA" = "Extraversion_(120k)",
      "FAMSAT" = "Family_Satisfaction_(160k)",
      "FRIENDSAT" = "Friend_Satisfaction_(160k)",
      #"FRIENDSAT1" = "Friend_Satisfaction_1_(160k)",
      #"FRIENDSAT2" = "Friend_Satisfaction_2_(160k)",
      "HAYFEVERb" = "Hayfever_(450k)",
      "HEIGHT" = "Height_(700k)",
      "LEFTOUT" = "Left Out_of Social_Activity_(510k)",
      "MENARCHE" = "Age_First_Menses_(330k)",
      "MIGRAINE" = "Migraine_(680k)",
      #"NEARSIGHTED" = "Nearsighted_(150k)",
      #"NEBWOM_UKBonly" = "Number_Ever Born_(women, UKB-only)_(240k)",
      "NEBwomen" = "Number_Ever Born_(women)_(240k, 470k)",
      "NEURO" = "Neuroticism_(480k)",
      "OPEN" = "Openness_(60k)",
      "RELIGATT" = "Religious_Attendance_(440k)",
      "RISK" = "Risk_(1430k)",
      "SELFHEALTH" = "Self-Rated_Health_(1200k)",
      "SWB" = "Subjective_Well-Being_(1020k)"
    ),
    phenotype_n = phenotype
  ) %>%
  filter(!(phenotype %in% unused_phenotypes))

########################################################
########### Save data for split-chart script ###########
########################################################

fwrite(data, paste0("output/prediction/HRSandWLS_phenotypes_r2_cleaned.txt"))


########################################################
################## Annotation function #################
########################################################

add.line <- function(expression){
  gsub('_', '\n', expression)
}


########################################################
################ Separate data by scale ################
########################################################

data_left <- data %>%
  filter(phenotype != "BMI_(760k)__" & phenotype != "EA_(1050k)" & phenotype != "Height_(700k)")

data_right <- data %>%
  filter(phenotype == "BMI_(760k)__" | phenotype == "EA_(1050k)" | phenotype == "Height_(700k)")


########################################################
###################### Plot chart ######################
########################################################

plot_left <- ggplot(data_left, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(hjust=0.65),
        axis.text = element_text(size = 10.5, colour="grey12"),
        legend.position = "none") +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data_left$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,15,2),
                     limits = c(0, 10.5),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2))) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

plot_right <- ggplot(data_right, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data_right$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,30,5),
                     limits = c(0, 31.5),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

grob <- arrangeGrob(
    plot_left,
    plot_right,

    nrow=1,
    widths=c(10, 2)
  )

grid.arrange(grob)

ggsave(
  file = paste0(
    "output/prediction/incr2_single_",
    Sys.Date(),
    ".png"
  ),
  plot = grob,
  width = 26,
  height = 8
)
