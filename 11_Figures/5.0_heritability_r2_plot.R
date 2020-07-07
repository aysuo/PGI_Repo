#----------------------------------------------------------------------------------#
# Charts r2 against heritability
# Date: 04/14/2020
# Author: Joel Becker

# Notes:
#
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "dplyr", "ggplot2", "stringr")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# source and set directory
source("directory_paths.R")
setwd(joel_PGS_dir)


########################################################
##################### Load r2 data #####################
########################################################

r2_path <- "output/prediction/HRSandWLS_phenotypes_r2_cleaned.txt"
r2_data <- fread(r2_path)


########################################################
################ Load heritability data ################
########################################################

unique(r2_data$phenotype)

h2_data <- data.frame(
  phenotype = unique(r2_data$phenotype),
  h2 = c(
    5.6,   # SWB
    4.39,  # risk
    4.71,  # relig (social)
    10.53, # openness
    11.22, # neuroticism
    NA,    # NEBWOM UKB only
    7.3,   # migraine
    48.99, # height
    11.68, # hayfever
    17.94, # extraversion
    9.1,   # ever smoker
    11.24, # educational attainment
    6.71,  # depression
    8.99,  # cigs per day
    22.53, # BMI
    5.4,   # asthma
    8.05,  # adventurousness
    8.19,  # self-rated health
    5.99,  # number ever born women
    5.85,  # left out social
    21.73, # intelligence
    19.53, # asthma ecz rhi
    19.36, # age first menses
    19.37, # age first birth
    21.92  # adhd
  ),
  h2_se = c(
    0.0020, # SWB
    0.0021, # risk
    0.0021, # relig (social)
    0.0090, # openness
    0.0042, # neuroticism
    NA,     # NEBWOM UKB only
    0.0031, # migraine
    0.0229, # height
    0.0071, # hayfever
    0.0096, # extraversion
    0.0030, # ever smoker
    0.0026, # educational attainment
    0.0025, # depression
    0.0085, # cigs per day
    0.0061, # BMI
    0.0047, # asthma
    0.0027, # adventurousness
    0.0019, # self-rated health
    0.0033, # number ever born women
    0.0022, # left out social
    0.0069, # intelligence
    0.0133, # asthma ecz rhi
    0.0071, # age first menses
    0.0073, # age first birth
    0.0127  # adhd
  ),
  pheno_category = c(
    "Well-being",     # SWB
    "Personality",    # risk
    "Personality",    # relig (social)
    "Personality",    # openness
    "Personality",    # neuroticism
    "Fertility",      # NEBWOM UKB only
    "Health",         # migraine
    "Anthropometric", # height
    "Health",         # hayfever
    "Personality",    # extraversion
    "Health",         # ever smoker
    "Cognition",      # educational attainment
    "Psychiatric",    # depression
    "Health",         # cigs per day
    "Anthropometric", # BMI
    "Health",         # asthma
    "Personality",    # adventurousness
    "Health",         # self-rated health
    "Fertility",      # number ever born women
    "Well-being",     # left out social
    "Cognition",      # intelligence
    "Health",         # asthma ecz rhi
    "Fertility",      # age first menses
    "Fertility",      # age first birth
    "Psychiatric"     # adhd
  )
)

# confirm length or phenotype list and h2 equal
# test would fail if phenotypes added/subtracted since I assigned h2
n_phenotype <- length(h2_data$phenotype[!is.na(h2_data$phenotype)])
n_h2 <- length(h2_data$h2[!is.na(h2_data$h2)])
n_h2_se <- length(h2_data$h2_se[!is.na(h2_data$h2_se)])
print(paste0("Number of phenotypes in h2 table: ", n_phenotype))
print(paste0("Number of h2 values in h2 table: ", n_h2))
print(paste0("Number of h2 se values in h2 table: ", n_h2_se))


########################################################
###################### Merge data ######################
########################################################

data <- r2_data %>%
  full_join(h2_data, by="phenotype")


########################################################
################# Plot chart functions #################
########################################################

r2_h2_plot <- function(data = data, cohort = c("HRS", "WLS"), error_bar_width = 0.5) {
  # creates plot of r2 against h2

  # build core of plot, depending on number of cohorts
  if (length(cohort) > 1) {
    # add deterministic jitter to horizontal axis data
    data <- data %>% mutate(h2 = case_when(
      Dataset == "HRS" ~ h2 + 0.25,
      Dataset == "WLS" ~ h2 - 0.25
    ))

    # build core plot
    plot <- ggplot(
        data,
        aes(fill = Dataset, colour = Dataset, x = h2, y = r2_inc)
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "longdash",
        colour = "#b4b4b4"
      ) +
      geom_errorbar(
        aes(ymin = r2_inc_lower, ymax = r2_inc_upper),
        width = error_bar_width
      ) +
      geom_errorbarh(
        aes(xmin = h2 - (1.96 * h2_se * 100), xmax = h2 + (1.96 * h2_se * 100)),
        height = 0.65 * error_bar_width
      ) +
      geom_point() +
      scale_color_manual(values = c("#1b9e77", "#7570b3"))
  } else if (length(cohort) == 1) {
    # filter data by cohort
    data_single_cohort <- filter(data, Dataset == cohort)
    # TODO: instead of creating separate plots, make this a facet plot if Dan wants?

    # build core plot
    colours <- c(
      "#1b9e77",
      "#d95f02",
      "#7570b3",
      "#e7298a",
      "#66a61e",
      "#e6ab02",
      "#a6761d"
    )
    plot <- ggplot(
        data_single_cohort,
        aes(colour = pheno_category, x = h2, y = r2_inc)
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "longdash",
        colour = "#b4b4b4"
      ) +
      geom_errorbar(
        aes(ymin = r2_inc_lower, ymax = r2_inc_upper),
        width = error_bar_width
      ) +
      geom_errorbarh(
        aes(xmin = h2 - (1.96 * h2_se * 100), xmax = h2 + (1.96 * h2_se * 100)),
        height = 0.65 * error_bar_width
      ) +
      geom_point() +
      scale_color_manual(values = colours) +
      labs(colour = "Phenotype\ncategories")
  }

  # build remainder of plot
  plot <- plot +
    scale_y_continuous(breaks = seq(0, 30, 5),
                       limits = c(0, max(data$r2_inc_upper) * 1.02),
                       expand = expand_scale(mult = c(0.01, 0), add = c(0, 0)),
                       labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = seq(0, 55, 5),
                       limits = c(0, max(data$h2) * 1.02),
                       expand = expand_scale(mult = c(0.01, 0), add = c(0, 1)),
                       labels = function(x) paste0(x, "%")) +
    theme_minimal() +
    ylab(expression(paste("Incremental ", italic(R)^2))) +
    xlab("Heritability")

  # return plot
  return(plot)
}

save_r2_h2_plot <- function(plot, width = 6, height = 5) {
  # saves plot, with name depending on plot input

  cohort <- deparse(substitute(plot))
  cohort <- str_split(cohort, "_plot", simplify = TRUE)[1,1]

  ggsave(
    file = paste0(
      "output/heritability/r2_h2_",
      cohort,
      "_",
      Sys.Date(),
      ".png"
    ),
    plot = plot,
    width = width,
    height = height
  )
}


########################################################
################# Plot and save charts #################
########################################################

HRS_plot <- r2_h2_plot(data = data, cohort = "HRS")
save_r2_h2_plot(plot = HRS_plot)

WLS_plot <- r2_h2_plot(data = data, cohort = "WLS")
save_r2_h2_plot(plot = WLS_plot)
