#----------------------------------------------------------------------------------#
# Compares prediction between public and single-trait scores graphically
# Date: 05/13/2020
# Author: Joel Becker

# Notes:
#  TODO: make differences more explicit and easier to handle by ggplot
#        by having *3* types of data; public, single positive, single
#        negative; then would show up fine in colours and legend
#  TODO: multiple fill scales, as here, separating cohorts and scores
#        https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
#----------------------------------------------------------------------------------#


########################################################
######################## Set-up ########################
########################################################

# load libraries
packages <- c("data.table", "tidyr", "dplyr", "ggplot2", "gridExtra")
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

# source and set directory
source("directory_paths.R")
setwd(joel_PGS_dir)


########################################################
####################### Load data ######################
########################################################

data_comparison <- fread(paste0(
  getwd(),
  "/output/prediction/public_single_comparison_r2.txt")
)


########################################################
##################### Wrangle data #####################
########################################################

# eventually want data to look like:
# phenotype | cohort | score_category | inc_r2

# but inc_r2 has to be manipulated for plot
# https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html

# inc_r2 for public scores:
# inc_r2 = abs(inc_r2_single - inc_r2_public)

# inc_r2 for single scores:
# if (inc_r2_single >= inc_r2_public) {
#   inc_r2 = inc_r2_single
# } else if (inc_r2_single < inc_r2_public) {
#   inc_r2 = inc_r2_single - inc_r2_public
# }

# for now, no SEs
barwidth <- 0.6
x_stretch_factor <- 1.5
HRS_shift_factor <- 0
WLS_shift_factor <- barwidth

# wrangle data for plot
select(data_comparison, phenotype, r2_inc_single, r2_inc_public)
data_plot <- data_comparison %>%
  mutate(
    single_greaterthan_public = case_when(
      r2_inc_single >= r2_inc_public ~ 1,
      r2_inc_single <  r2_inc_public ~ 0
    ),
    #incr2_single_positive_wrangled = case_when(
    #  single_greaterthan_public == 1 ~ r2_inc_single - r2_inc_public,
    #  TRUE ~ 0
    #),
    #incr2_single_negative_wrangled = case_when(
    #  single_greaterthan_public == 0 ~ -1 * (r2_inc_single - r2_inc_public),
    #  TRUE ~ 0
    #),
    incr2_single_positive_wrangled_HRS = case_when(
      single_greaterthan_public == 1 & Dataset == "HRS" ~ r2_inc_single - r2_inc_public,
      TRUE ~ 0
    ),
    incr2_single_negative_wrangled_HRS = case_when(
      single_greaterthan_public == 0 & Dataset == "HRS" ~ -1 * (r2_inc_single - r2_inc_public),
      TRUE ~ 0
    ),
    incr2_single_positive_wrangled_WLS = case_when(
      single_greaterthan_public == 1 & Dataset == "WLS" ~ r2_inc_single - r2_inc_public,
      TRUE ~ 0
    ),
    incr2_single_negative_wrangled_WLS = case_when(
      single_greaterthan_public == 0 & Dataset == "WLS" ~ -1 * (r2_inc_single - r2_inc_public),
      TRUE ~ 0
    ),
    incr2_public_HRS_wrangled = case_when(
      Dataset == "HRS" ~ r2_inc_public,
      TRUE ~ 0
    ),
    incr2_public_WLS_wrangled = case_when(
      Dataset == "WLS" ~ r2_inc_public,
      TRUE ~ 0
    )
  )

# reshape data
data_plot <- data_plot %>%
  select(
    phenotype,
    cohort = Dataset,
    #incr2_single = incr2_single_wrangled,
    incr2_publicHRS = incr2_public_HRS_wrangled,
    incr2_publicWLS = incr2_public_WLS_wrangled,
    #incr2_singlepositive = incr2_single_positive_wrangled,
    #incr2_singlenegative = incr2_single_negative_wrangled,
    incr2_singlepositiveHRS = incr2_single_positive_wrangled_HRS,
    incr2_singlenegativeHRS = incr2_single_negative_wrangled_HRS,
    incr2_singlepositiveWLS = incr2_single_positive_wrangled_WLS,
    incr2_singlenegativeWLS = incr2_single_negative_wrangled_WLS,
    single_greaterthan_public
  ) %>%
  gather(
    key = "score_category",
    value = "inc_r2",
    incr2_publicHRS,
    incr2_publicWLS,
    #incr2_singlepositive,
    #incr2_singlenegative,
    incr2_singlepositiveHRS,
    incr2_singlenegativeHRS,
    incr2_singlepositiveWLS,
    incr2_singlenegativeWLS
  ) %>%
  separate("score_category", c("incr2_separated", "score_category")) %>%
  select(-"incr2_separated") %>%
  group_by(phenotype) %>%
  mutate(
    n_cohorts = length(unique(cohort)),
    cohorts_available = case_when(
      n_cohorts == 2 ~ "both",
      n_cohorts == 1 & cohort == "HRS" ~ "HRS",
      n_cohorts == 1 & cohort == "WLS" ~ "WLS"
    )
  ) %>%
  arrange(phenotype, cohort, score_category)

# phenotype labels
phenotype_labels <- data.frame(
  phenotype = unique(data_plot[, "phenotype"]),
  phenotype_label = c(
    "ADHD",
    "Adventurousness",
    "Age_First_Birth",
    "Asthma_Eczema_Rhinitis",
    "Asthma",
    "BMI___",
    "Cognitive_Performance",
    "Cigarettes_per Day",
    "Depression",
    #"Educational_Attainment__",
    "Ever_Smoker",
    "Extraversion",
    "Family_Satisfaction_",
    "Friend_Satisfaction_",
    "Friend_Satisfaction_1",
    "Friend_Satisfaction_2",
    "Hayfever",
    "Height___",
    "Left_Out of_Social_Activity",
    "Menarche",
    "Migraine",
    "Number_ever born_(women)",
    "Neuroticism",
    "Openness",
    "Religious_Attendance",
    "Risk",
    "Self-_Rated_Health",
    "Subjective_Well-_Being"
  )
)

########################################################
############## Organise data labels, left ##############
########################################################

# details on conditional bar colouring here
# https://stackoverflow.com/questions/36613033/conditionally-colouring-bars-in-ggplot2

# stack and dodge is hard but doable, see here
# https://stackoverflow.com/questions/12715635/ggplot2-bar-plot-with-both-stack-and-dodge
# or better solution here
# https://community.rstudio.com/t/ggplot-position-dodge-with-position-stack/16425/7

data_left <- data_plot %>%
  filter(
    phenotype != "BMI" &
    phenotype != "EA" &
    phenotype != "HEIGHT"
  )

# index phenotypes
phenotype_x_left <- data.frame(
    phenotype = data_left$phenotype,
    phenotype_index = data_left %>% group_indices(phenotype),
    n_cohorts = data_left$n_cohorts,
    cohorts_available = data_left$cohorts_available
  ) %>%
  mutate(
    phenotype_index = phenotype_index * x_stretch_factor,
    phenotype_x = phenotype_index - (min(phenotype_index) - 0.1)
  ) %>%
  distinct(.)

for (row in 1:nrow(phenotype_x_left)) {
  # shift all rows left if only one cohort available
  if (phenotype_x_left[row, "cohorts_available"] == "HRS") {
    for (inner_row in row:nrow(phenotype_x_left)) {
      phenotype_x_left[inner_row, "phenotype_x"] <- phenotype_x_left[inner_row, "phenotype_x"] - HRS_shift_factor
    }
  } else if (phenotype_x_left[row, "cohorts_available"] == "WLS") {
    for (inner_row in row:nrow(phenotype_x_left)) {
      phenotype_x_left[inner_row, "phenotype_x"] <- phenotype_x_left[inner_row, "phenotype_x"] - WLS_shift_factor
    }
  }
}

for (row in 1:(nrow(phenotype_x_left) - 1)) {
  # shift all rows left disproportionately if only one cohort available
  if (phenotype_x_left[row, "cohorts_available"] == "HRS") {
    if (phenotype_x_left[(row + 1), "cohorts_available"] == "WLS") {
      for (inner_row in (row + 1):nrow(phenotype_x_left)) {
        phenotype_x_left[inner_row, "phenotype_x"] <- phenotype_x_left[inner_row, "phenotype_x"] - WLS_shift_factor
      }
    } else if (phenotype_x_left[(row + 1), "cohorts_available"] == "both") {
      for (inner_row in (row + 1):nrow(phenotype_x_left)) {
        phenotype_x_left[inner_row, "phenotype_x"] <- phenotype_x_left[inner_row, "phenotype_x"] - WLS_shift_factor
      }
    }
  }
}

for (row in 1:nrow(phenotype_x_left)) {
  # get left-/right-most x-values for illustrative purposes
  if (phenotype_x_left[row, "n_cohorts"] == 1) {
    phenotype_x_left[row, "phenotype_x_lower"] <- phenotype_x_left[row, "phenotype_x"] - (barwidth / 2)
      phenotype_x_left[row, "phenotype_x_upper"] <- phenotype_x_left[row, "phenotype_x"] + (barwidth / 2)
  } else if (phenotype_x_left[row, "n_cohorts"] == 2) {
    phenotype_x_left[row, "phenotype_x_lower"] <- phenotype_x_left[row, "phenotype_x"] - barwidth
      phenotype_x_left[row, "phenotype_x_upper"] <- phenotype_x_left[row, "phenotype_x"] + barwidth
  }

  # shift labels slightly if only one cohort available
  if (phenotype_x_left[row, "cohorts_available"] == "HRS") {
    phenotype_x_left[row, "phenotype_label_x"] <- phenotype_x_left[row, "phenotype_x"] - (barwidth / 2)
  } else if (phenotype_x_left[row, "cohorts_available"] == "WLS") {
    phenotype_x_left[row, "phenotype_label_x"] <- phenotype_x_left[row, "phenotype_x"] + (barwidth / 2)
  } else if (phenotype_x_left[row, "cohorts_available"] == "both") {
    phenotype_x_left[row, "phenotype_label_x"] <- phenotype_x_left[row, "phenotype_x"]
  }
}

data_left <- data_left %>%
  full_join(
    phenotype_x_left,
    by = c("phenotype", "n_cohorts")
  )

#select(data_plot, phenotype, cohort, score_category, inc_r2, phenotype_x_lower, phenotype_x_upper)

phenotype_x_left_labels <- data_left %>%
  select(phenotype, phenotype_x, phenotype_label_x) %>%
  distinct(.) %>%
  left_join(phenotype_labels, by = "phenotype")

data_left_HRS <- data_left %>%
  filter(cohort == "HRS") %>%
  group_by(phenotype) %>%
  arrange(score_category)

data_left_WLS <- data_left %>%
  filter(cohort == "WLS") %>%
  group_by(phenotype) %>%
  arrange(score_category)


########################################################
############## Organise data labels, right #############
########################################################

data_right <- data_plot %>%
  filter(
    phenotype == "BMI" |
    phenotype == "EA" |
    phenotype == "HEIGHT"
  )

# index phenotypes
phenotype_x_right <- data.frame(
    phenotype = data_right$phenotype,
    phenotype_index = data_right %>% group_indices(phenotype),
    n_cohorts = data_right$n_cohorts,
    cohorts_available = data_right$cohorts_available
  ) %>%
  mutate(
    phenotype_index = phenotype_index * x_stretch_factor,
    phenotype_x = phenotype_index - (min(phenotype_index) - 0.1)
  ) %>%
  distinct(.)

for (row in 1:nrow(phenotype_x_right)) {
  # shift all rows left if only one cohort available
  if (phenotype_x_right[row, "cohorts_available"] == "HRS") {
    for (inner_row in row:nrow(phenotype_x_right)) {
      phenotype_x_right[inner_row, "phenotype_x"] <- phenotype_x_right[inner_row, "phenotype_x"] - HRS_shift_factor
    }
  } else if (phenotype_x_right[row, "cohorts_available"] == "WLS") {
    for (inner_row in row:nrow(phenotype_x_right)) {
      phenotype_x_right[inner_row, "phenotype_x"] <- phenotype_x_right[inner_row, "phenotype_x"] - WLS_shift_factor
    }
  }
}

for (row in 1:(nrow(phenotype_x_right) - 1)) {
  # shift all rows left disproportionately if only one cohort available
  if (phenotype_x_right[row, "cohorts_available"] == "HRS") {
    if (phenotype_x_right[(row + 1), "cohorts_available"] == "WLS") {
      for (inner_row in (row + 1):nrow(phenotype_x_right)) {
        phenotype_x_right[inner_row, "phenotype_x"] <- phenotype_x_right[inner_row, "phenotype_x"] - WLS_shift_factor
      }
    } else if (phenotype_x_right[(row + 1), "cohorts_available"] == "both") {
      for (inner_row in (row + 1):nrow(phenotype_x_right)) {
        phenotype_x_right[inner_row, "phenotype_x"] <- phenotype_x_right[inner_row, "phenotype_x"] - WLS_shift_factor
      }
    }
  }
}

for (row in 1:nrow(phenotype_x_right)) {
  # get left-/right-most x-values for illustrative purposes
  if (phenotype_x_right[row, "n_cohorts"] == 1) {
    phenotype_x_right[row, "phenotype_x_lower"] <- phenotype_x_right[row, "phenotype_x"] - (barwidth / 2)
      phenotype_x_right[row, "phenotype_x_upper"] <- phenotype_x_right[row, "phenotype_x"] + (barwidth / 2)
  } else if (phenotype_x_right[row, "n_cohorts"] == 2) {
    phenotype_x_right[row, "phenotype_x_lower"] <- phenotype_x_right[row, "phenotype_x"] - barwidth
      phenotype_x_right[row, "phenotype_x_upper"] <- phenotype_x_right[row, "phenotype_x"] + barwidth
  }

  #  shift labels slightly if only one cohort available
  if (phenotype_x_right[row, "cohorts_available"] == "HRS") {
    phenotype_x_right[row, "phenotype_label_x"] <- phenotype_x_right[row, "phenotype_x"] - (barwidth / 2)
  } else if (phenotype_x_left[row, "cohorts_available"] == "WLS") {
    phenotype_x_right[row, "phenotype_label_x"] <- phenotype_x_right[row, "phenotype_x"] + (barwidth / 2)
  } else if (phenotype_x_left[row, "cohorts_available"] == "both") {
    phenotype_x_right[row, "phenotype_label_x"] <- phenotype_x_right[row, "phenotype_x"]
  }
}

data_right <- data_right %>%
  full_join(
    phenotype_x_right,
    by = c("phenotype", "n_cohorts")
  )

#select(data_plot, phenotype, cohort, score_category, inc_r2, phenotype_x_lower, phenotype_x_upper)

phenotype_x_right_labels <- data_right %>%
  select(phenotype, phenotype_x, phenotype_label_x) %>%
  distinct(.) %>%
  left_join(phenotype_labels, by = "phenotype")

data_right_HRS <- data_right %>%
  filter(cohort == "HRS") %>%
  group_by(phenotype) %>%
  arrange(score_category)

data_right_WLS <- data_right %>%
  filter(cohort == "WLS") %>%
  group_by(phenotype) %>%
  arrange(score_category)


########################################################
####################### Plot data ######################
########################################################

add.line <- function(expression){
  gsub('_', '\n', expression)
}

plot_left <- ggplot() +
    geom_bar(data = data_left_HRS,
             mapping = aes(
               x = phenotype_x - (barwidth/2),
               y = inc_r2,
               fill = factor(
                 score_category,
                 levels=c("singlepositiveHRS", "singlenegativeHRS", "singlepositiveWLS", "singlenegativeWLS", "publicHRS", "publicWLS")
               )
             ),
             colour = "black",
             stat="identity",
             position='stack',
             width = barwidth) +
    #geom_text(data = month_one,
    #          aes(x = id, y = pos, label = count )) +
    geom_bar(data = data_left_WLS,
             mapping = aes(
               x = phenotype_x + (barwidth/2),
               y = inc_r2,
               fill = factor(
                 score_category,
                 levels=c("singlepositiveHRS", "singlenegativeHRS", "singlepositiveWLS", "singlenegativeWLS", "publicHRS", "publicWLS")
               )
             ),
             colour = "black",
             stat="identity",
             position='stack' ,
             width = barwidth) +
    #geom_text(data = month_two,
    #          aes(x = id + barwidth + 0.01, y = pos, label = count )) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 17, colour="grey4"),
      axis.title.x = element_text(hjust=0.65),
      axis.text = element_text(size = 10.5, colour="grey12"),
      legend.position = "none"
    ) +
    scale_y_continuous(
      breaks = seq(0, 8, 2),
      limits = c(0, 8),
      expand = c(0, 0),
      labels = function(x) paste0(x, "%")
    ) +
    scale_x_continuous(
      breaks = c(phenotype_x_left_labels$phenotype_label_x),
      expand = c(0, 0),
      label = add.line(c(as.character(phenotype_x_left_labels$phenotype_label)))
    ) +
    xlab("Phenotype") +
    ylab(expression(paste("Incremental ", italic(R)^2))) +
    scale_fill_manual(
      values = c(
        "singlepositiveHRS" = "#006d2c",
        "singlenegativeHRS" = "#bae4b3",
        "singlepositiveWLS" = "#54278f",
        "singlenegativeWLS" = "#cbc9e2",
        "publicHRS" = "#31a354",
        "publicWLS" = "#756bb1"
      ),
      name = "",
      labels = c(
        "Improvement (HRS, positive)",
        "Improvement (HRS, negative)",
        "Improvement (WLS, positive)",
        "Improvement (WLS, negative)",
        "Public, HRS",
        "Public, WLS"
      )
    )

plot_right <- ggplot() +
    geom_bar(data = data_right_HRS,
             mapping = aes(
               x = phenotype_x - (barwidth/2),
               y = inc_r2,
               fill = factor(
                 score_category,
                 levels=c("singlepositiveHRS", "singlenegativeHRS", "singlepositiveWLS", "singlenegativeWLS", "publicHRS", "publicWLS")
               )
             ),
             colour = "black",
             stat="identity",
             position='stack',
             width = barwidth) +
    geom_bar(data = data_right_WLS,
             mapping = aes(
               x = phenotype_x + (barwidth/2),
               y = inc_r2,
               fill = factor(
                 score_category,
                 levels=c("singlepositiveHRS", "singlenegativeHRS", "singlepositiveWLS", "singlenegativeWLS", "publicHRS", "publicWLS")
               )#,
               #colour = c("#1b9e77", "#7570b3")
             ),
             colour = "black",
             stat="identity",
             position='stack' ,
             width = barwidth) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 17, colour="grey4"),
      #axis.title.x = element_text(),
      axis.text = element_text(size = 10.5, colour="grey12")
    ) +
    scale_y_continuous(
      breaks = seq(0,30,5),
      limits = c(0, 30),
      expand = c(0, 0),
      labels = function(x) paste0(x, "%")
    ) +
    scale_x_continuous(
      breaks = c(phenotype_x_right_labels$phenotype_label_x),
      expand = c(0, 0),
      label = add.line(c(as.character(phenotype_x_right_labels$phenotype_label)))
    ) +
    xlab("") +
    ylab("") +
    scale_fill_manual(
      values = c(
        "singlepositiveHRS" = "#006d2c",
        "singlenegativeHRS" = "#bae4b3",
        "singlepositiveWLS" = "#54278f",
        "singlenegativeWLS" = "#cbc9e2",
        "publicHRS" = "#31a354",
        "publicWLS" = "#756bb1"
      ),
      name = "",
      labels = c(
        "Improvement (HRS, positive)",
        "Improvement (HRS, negative)",
        "Improvement (WLS, positive)",
        "Improvement (WLS, negative)",
        "Public, HRS",
        "Public, WLS"
      )
    )

grob <- arrangeGrob(
    plot_left,
    plot_right,

    nrow=1,
    widths=c(10, 2.3)
  )

#grid.arrange(grob)

ggsave(
  file = paste0(
    "output/prediction/incr2_single_public_comparison_",
    Sys.Date(),
    ".png"
  ),
  plot = grob,
  width = 26,
  height = 10
)








# stack and dodge code
# separate the by the variable which we are dodging by so
# we have two data frames impact and control
data_plot <- data_plot[1:18, ]
data_left_HRS <- data_plot %>% filter(cohort == "HRS")
#impact <- df %>% filter(treatment == "Impact") %>%
#  mutate(pos = sum(total, na.rm=T))

data_left_WLS <- data_plot %>% filter(cohort == "WLS")
#control <- df %>% filter(treatment == "Control") %>%
#  mutate(pos = sum(total, na.rm=T))

# calculate the position for the annotation element
#impact_an <- impact %>% group_by(year) %>%
#  summarise(
#    pos = sum(total) + 12
#    , treatment = first(treatment)
#  )

#control_an <- control %>% group_by(year) %>%
#  summarise(
#    pos = sum(total) + 12
#    , treatment = first(treatment)
#  )

# define the width of the bars, we need this set so that
# we can use it to position the second layer geom_bar
barwidth = 0.8

ggplot() +
  geom_bar(
    data = data_left_HRS
    , aes(x = phenotype_x - ((barwidth + 0.1) / 2), y = inc_r2, fill = score_category)
    , position = "stack"
    , stat = "identity"
    , width = barwidth
  ) +
  #geom_bar(
  #  data = impact
  #  , aes(x = year, y = total, fill = type)
  #  , position = "stack"
  #  , stat = "identity"
  #  , width = barwidth
  #) +
  #annotate(
  #  "text"
  #  , x = impact_an$year
  #  ,y = impact_an$pos
  #  , angle = 90
  #  , label = impact_an$treatment
  #) +
  geom_bar(
    data = data_left_WLS
    # here we are offsetting the position of the second layer bar
    # by adding the barwidth plus 0.1 to push it to the right
    , aes(x = phenotype_x + ((barwidth + 0.1) / 2), y = inc_r2, fill = score_category)
    , position = "stack"
    , stat = "identity"
    , width = barwidth
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 17, colour="grey4"),
    axis.title.x = element_text(),
    axis.text = element_text(size = 10.5, colour="grey12"),
    legend.position = "none"
  ) +
  scale_y_continuous(
    #breaks = seq(0,15,2),
    #limits = c(0, 12.6),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2)))
  #geom_bar(
  #  data = control
  #  # here we are offsetting the position of the second layer bar
  #  # by adding the barwidth plus 0.1 to push it to the right
  #  , aes(x = year + barwidth + 0.1, y = total, fill = type)
  #  , position = "stack"
  #  , stat = "identity"
  #  , width = barwidth
  #) +
  #annotate(
  #  "text"
  #  , x = control_an$year + (barwidth * 1) + 0.1
  #  ,y = control_an$pos
  #  , angle = 90
  #  , label = control_an$treatment
  #) +
  #scale_x_discrete(limits = c(2010, 2011))

# my plot code
plot <- ggplot(
    data_plot,
    aes(fill = score_category, y = inc_r2, x = phenotype)
  ) +
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 17, colour="grey4"),
    axis.title.x = element_text(),
    axis.text = element_text(size = 10.5, colour="grey12"),
    legend.position = "none"
  ) +
  geom_col(
    colour = "black",
    width = 0.9,
    position=position_dodge2(preserve = "single", padding = 0)
  ) +
  scale_x_discrete(labels = add.line(sort(unique(data_left$phenotype_n)))) +
  scale_y_continuous(
    #breaks = seq(0,15,2),
    #limits = c(0, 12.6),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  xlab("Phenotype") +
  ylab(expression(paste("Incremental ", italic(R)^2)))
