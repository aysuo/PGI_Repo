#----------------------------------------------------------------------------------#
# Charts incremental r2 of phenotypes, split in two for ease of reading
# Date: 03/02/2020
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

path <- "output/prediction/HRSandWLS_phenotypes_r2_cleaned.txt"
data <- fread(path) %>%
  arrange(phenotype)


########################################################
################## Annotation function #################
########################################################

add.line <- function(expression){
  gsub('_', '\n', expression)
}


########################################################
############# Separate data by plots/scale #############
########################################################

# separate right-most phenotypes
data_plot2b <- data %>%
  filter(phenotype == "BMI_(760k)__" |
         phenotype == "EA_(1050k)" |
         phenotype == "Height_(700k)")

data_leftplot <- data %>%
  anti_join(data_plot2b)

# separate out left plot vs right_a plot phenotypes
unique_left_phenotypes <- length(unique(data_leftplot$phenotype))
unique_right_phenotypes <- length(unique(data_plot2b$phenotype))
n_phenotypes_plot1 <- (unique_left_phenotypes + unique_right_phenotypes + 2) %/% 2

phenotypes_plot1 <- unique(data_leftplot$phenotype)[1:n_phenotypes_plot1]
phenotypes_plot2a <- data_leftplot %>%
  filter(!(phenotype %in% phenotypes_plot1))
phenotypes_plot2a <- unique(phenotypes_plot2a$phenotype)

data_plot1 <- data %>%
  filter(phenotype %in% phenotypes_plot1)

data_plot2a <- data %>%
  filter(phenotype %in% phenotypes_plot2a)


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
                     limits = c(0, 10.5),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2))) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

ggsave(
  file = paste0(
    "output/prediction/incr2_single_split_left_",
    Sys.Date(),
    ".png"
  ),
  plot = plot1,
  width = 14,
  height = 7
)

plot2a <- ggplot(data_plot2a, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.title.x = element_text(hjust=0.7),
        axis.text = element_text(size = 10.5, colour="grey12"),
        legend.position = "none") +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data_plot2a$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,15,2),
                     limits = c(0, 10.5),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2))) +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

plot2b <- ggplot(data_plot2b, aes(fill=Dataset, x = phenotype, y = r2_inc)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 17, colour="grey4"),
        axis.text = element_text(size = 10.5, colour="grey12")) +
  geom_col(colour = "black", width = 0.9, position=position_dodge2(preserve = "single", padding = 0)) +
  geom_errorbar(aes(ymin = r2_inc_lower, ymax = r2_inc_upper), width = .1, position=position_dodge(width = 0.9)) +
  scale_x_discrete(labels = add.line(sort(unique(data_plot2b$phenotype_n)))) +
  scale_y_continuous(breaks = seq(0,30,5),
                     limits = c(0, 31.5),
                     expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))

grob <- arrangeGrob(
    plot2a,
    plot2b,

    nrow=1,
    widths=c(5.4, 2.6)
  )

grid.arrange(grob)

ggsave(
  file = paste0(
    "output/prediction/incr2_single_split_right_",
    Sys.Date(),
    ".png"
  ),
  plot = grob,
  width = 14,
  height = 7
)
