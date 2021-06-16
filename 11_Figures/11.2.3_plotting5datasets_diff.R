#----------------------------------------------------------------------------------#
# Charts differences in single trait and multi trait results
# Date: 04/16/2019
# Author: Joel Becker, Hariharan Jayashankar
#----------------------------------------------------------------------------------#

# load libraries
packages <- c("data.table", "dplyr", "ggplot2", 
              "gridExtra", "grid", "gtable", "rlang",
              "stringr", "viridis", "cowplot", "latex2exp") #, "pBrackets"
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = TRUE)

args=commandArgs(trailingOnly=TRUE)
mainDir=args[1]

datapath <- paste0(mainDir,"/derived_data/11_Figures/Prediction/input/")
outpath <- paste0(mainDir,"/derived_data/11_Figures/Prediction/output/")

########################################################
####################### Load data ######################
########################################################

data <- readRDS(paste0(datapath, "phenotypes_r2_cleaned.R"))
phenotype_labels <- fread(paste0(datapath, "phenotype_labels.txt"))

# Preparing for Plotting ----------------------------------------------------------------

########################################################
################# Annotation functions #################
########################################################

add.line <- function(expression){
    
    # Adds \n (line breaks) to x labels in graphs
    
    gsub('_', '\n', expression)
}

make.x.labels <- function(phenotype_name, single_or_multi_score){
    
    # replaces variable names with their labels
    
    # init empty vector
    pheno_label_vec <- c()
    
    # loop over vector of phenotype_names
    i = 1
    for (name in phenotype_name) {
        
        pheno_label <- phenotype_labels[phenotype == name, phenotype_label]       
        pheno_label <- add.line(pheno_label)
        # add to vector
        pheno_label_vec[[i]] <- pheno_label
        
        # iterate i
        i <- i +1
        
    }
    
    return(pheno_label_vec)
    
}


# Make variable which is equal to r2_inc_single if r2_diff_single_public is missing
data <- data %>%
    mutate(r2_diff_single_public2 = ifelse(complete.cases(r2_diff_single_public),
                                           r2_diff_single_public,
                                           r2_inc_single),
           r2_diff_multi_public2 = ifelse(complete.cases(r2_diff_multi_public),
                                          r2_diff_multi_public,
                                          r2_inc_multi))


# Make new group based on whether or not r2_diff_upper and lower are mising
# for different shaped for missing and non missing
data <- data %>%
    mutate(data_group2 = case_when(complete.cases(r2_diff_single_public) ~ "Available",
                                   TRUE ~ "Unavailable"))


data <- data %>%
    mutate(r2_diff_single_public_lower2 = ifelse(complete.cases(r2_diff_single_public_lower),
                                                 r2_diff_single_public_lower,
                                                 r2_inc_single_lower),
           r2_diff_single_public_upper2 = ifelse(complete.cases(r2_diff_single_public_upper),
                                                 r2_diff_single_public_upper,
                                                 r2_inc_single_upper))


data <- data %>%
    mutate(r2_diff_multi_public_lower2 = ifelse(complete.cases(r2_diff_multi_public_lower),
                                                r2_diff_multi_public_lower,
                                                r2_inc_multi_lower),
           r2_diff_multi_public_upper2 = ifelse(complete.cases(r2_diff_multi_public_upper),
                                                r2_diff_multi_public_upper,
                                                r2_inc_multi_upper))



########################################################
###################### Plot chart ######################
########################################################

facet_labels <- as_labeller(c(`Anthropometric` = "Anthropometric", 
                              `Cognition` = "Cognition and Education", 
                              `Fertility and Sexual Development` = "Fertility and Sexual Development",
                              `Health and Health Behaviors` = "Health and Health Behaviors",
                              `Personality and Well-Being` = "Personality and Well-Being"))

line <- linesGrob(unit(c(.05, .95), "npc"), unit(1, "npc"),
                  gp = gpar(col = 'darkgrey', lwd = 4))



trait_n <- "single"

# y axis var names
yvar_bar <- paste0("r2_inc_", trait_n)
yvar_bar_lb <- paste0("r2_inc_", trait_n, "_lower")
yvar_bar_ub <- paste0("r2_inc_", trait_n, "_upper")

yvar_diff <- paste0("r2_diff_", trait_n, "_public2")
yvar_diff_lb <- paste0("r2_diff_", trait_n, "_public_lower2")
yvar_diff_ub <- paste0("r2_diff_", trait_n, "_public_upper2")

groups <- data %>% pull(pheno_group) %>% unique %>% sort()

# specifying colors
colour_scale <- c("HRS" = "#1b9e77", "WLS" = "#7570b3", 
                  "UKB" = "#d95f02", "Dunedin" = "#d83951",
                  "E-Risk" = "#dbc020")

color_pal <- c("chartreuse3", "blue", "purple", "gold", "brown")

plot_212 <- data %>%
    filter(complete.cases(r2_diff_multi_single)) %>%
    filter(pheno_group %in% groups[1:3]) %>%
    filter(phenotype != "BMI" & phenotype != "HEIGHT") %>%
    ggplot(aes(x = phenotype_n, y = r2_diff_multi_single, group = Dataset)) +
    facet_grid( .~pheno_group ,  space = "free",  scales = "free", switch="x",
                labeller = facet_labels) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          axis.title = element_text(size = 17, colour="grey4"),
          axis.text.x = element_text(size = 10.5, colour="grey12"),
          axis.text.y = element_text(size = 12, 
                                     colour="grey12"),
          strip.placement = "outside",
          strip.text.x = element_text(size = 13, face = "bold"),
          panel.spacing = unit(0, "lines"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.position=c(.85, 0.78),
          legend.background = element_rect(colour = 'black',
                                           fill = 'white',
                                           linetype='solid'),
          legend.box = "horizontal")  +
    geom_errorbar(aes(ymin = r2_diff_multi_single_lower,
                      ymax = r2_diff_multi_single_upper,
                      color = Dataset),
                  width = .5,
                  size = 0.5,
                  position=position_dodge2(preserve = "single",
                                           width = 0.5)) +
    
    geom_point(aes(shape = Dataset,
                   color = Dataset),
               fill = "white",
               position=position_dodge2(width = 0.5, preserve = "single"),
               size = 3) +
    geom_point(aes(shape = Dataset, 
                   fill=Dataset,
                   color = Dataset), 
               position=position_dodge2(width = 0.5, preserve = "single"), 
               size = 3) +
    
    scale_x_discrete(labels = setNames(make.x.labels(data$phenotype_n, 
                                                     "sample_single"), 
                                       data$phenotype_n)) +
    scale_y_continuous(expand = c(0, 0),
                       labels = function(x) paste0(x, "%")) +
    geom_hline(yintercept = 0, color = "black") +
    xlab("") +
    ylab(TeX("$\\Delta$ Incremental $\\textit{R}^2$")) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    scale_shape_manual(values = 21:25) +
    coord_cartesian(clip = "off") +
    annotation_custom(line, ymin = -2.1, ymax = -2.1) +
    guides(color = guide_legend(nrow=3,byrow=TRUE),
           fill = guide_legend(nrow=3,byrow=TRUE),
           alpha = FALSE)

plot_22 <- data %>%
    filter(complete.cases(r2_diff_multi_single)) %>%
    filter(pheno_group == groups[4]) %>%
    ggplot(aes(x = phenotype_n, y = r2_diff_multi_single, group = Dataset)) +
    facet_grid( .~pheno_group,  space = "free",  scales = "free", switch="x",
                labeller = facet_labels) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          axis.title = element_text(size = 17, colour="grey4"),
          axis.text.x = element_text(size = 10.5, colour="grey12"),
          axis.text.y = element_text(size = 12, 
                                     colour="grey12"),
          strip.placement = "outside",
          strip.text.x = element_text(size = 13, face = "bold"),
          panel.spacing = unit(0, "lines"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position="none")  +
    geom_errorbar(aes(ymin = r2_diff_multi_single_lower, 
                      ymax = r2_diff_multi_single_upper,
                      color = Dataset), 
                  width = .5, 
                  size = 0.5,
                  position=position_dodge2(preserve = "single",
                                           width = 0.5)) +
    
    geom_point(aes(shape = Dataset,
                   color = Dataset),
               fill = "white",
               position=position_dodge2(width = 0.5, preserve = "single"),
               size = 3) +
    geom_point(aes(shape = Dataset, 
                   fill=Dataset,
                   color = Dataset), 
               position=position_dodge2(width = 0.5, preserve = "single"), 
               size = 3) +
    
    scale_x_discrete(labels = setNames(make.x.labels(data$phenotype_n, 
                                                     "sample_single"), 
                                       data$phenotype_n)) +
    scale_y_continuous(expand = c(0, 0),
                       labels = function(x) paste0(x, "%")) +
    geom_hline(yintercept = 0, color = "black") +
    xlab("") +
    ylab(TeX("$\\Delta$ Incremental $\\textit{R}^2$")) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    scale_shape_manual(values = 21:25) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_cartesian(clip = "off") +
    labs(shape = "Public Score") +
    annotation_custom(line, ymin = -2.1, ymax = -2.1) +
    guides(shape = FALSE)

plot_23 <- data %>%
    filter(complete.cases(r2_diff_multi_single)) %>%
    filter(pheno_group == groups[5]) %>%
    ggplot(aes(x = phenotype_n, y = r2_diff_multi_single, group = Dataset)) +
    facet_grid( .~pheno_group,  space = "free",  scales = "free", switch="x",
                labeller = facet_labels) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          axis.title = element_text(size = 17, colour="grey4"),
          axis.text.x = element_text(size = 10.5, colour="grey12"),
          axis.text.y = element_text(size = 12, 
                                     colour="grey12"),
          strip.placement = "outside",
          strip.text.x = element_text(size = 13, face = "bold"),
          panel.spacing = unit(0, "lines"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position="none")  +
    geom_errorbar(aes(ymin = r2_diff_multi_single_lower,
                      ymax = r2_diff_multi_single_upper,
                      color = Dataset),
                  width = .5,
                  size = 0.5,
                  position=position_dodge2(preserve = "single",
                                           width = 0.5)) +
    
    
    
    geom_point(aes(shape = Dataset,
                   color = Dataset),
               fill = "white",
               position=position_dodge2(width = 0.5, preserve = "single"),
               size = 3) +
    geom_point(aes(shape = Dataset, 
                   fill=Dataset,
                   color = Dataset), 
               position=position_dodge2(width = 0.5, preserve = "single"), 
               size = 3) +
    scale_x_discrete(labels = setNames(make.x.labels(data$phenotype_n, 
                                                     "sample_single"), 
                                       data$phenotype_n)) +
    scale_y_continuous(expand = c(0, 0),
                       labels = function(x) paste0(x, "%")) +
    geom_hline(yintercept = 0, color = "black") +
    xlab("") +
    ylab(TeX("$\\Delta$ Incremental $\\textit{R}^2$")) +
    scale_fill_manual(values = color_pal) +
    scale_color_manual(values = color_pal) +
    scale_shape_manual(values = 21:25) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    coord_cartesian(clip = "off") +
    labs(shape = "Public Score") +
    annotation_custom(line, ymin = -4.3, ymax = -4.3) +
    guides(shape = FALSE)


# arranging plots
g12 <- ggplotGrob(plot_212)
g12$grobs[[7]]$layout$clip <- "off"
g1 <- g12
g2 <- ggplotGrob(plot_22)
g3 <- ggplotGrob(plot_23)
g1$grobs[[7]]$layout$clip <- "off"
g2$grobs[[7]]$layout$clip <- "off"
g3$grobs[[7]]$layout$clip <- "off"
panelA <- grid.arrange(g1, g2, g3)


# outputting
ggsave(
    file = paste0(
        outpath,
        "diff_",
        "multi_single",
        ".png"
    ),
    plot = panelA,
    width = 16,
    height = 12
)

ggsave(
    file = paste0(
        outpath,
        "diff_",
        "multi_single",
        ".pdf"
    ),
    plot = panelA,
    width = 16,
    height = 12
)
