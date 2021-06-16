# ========================
# Drawing a Forest Plot of
# PGI Estimates
# ========================

library(tidyverse)
library(grid)
library(scales)
library(gridExtra)
library(ggpubr)
library(ggforestplot)
library(ggforce)
theme_set(theme_minimal())

# Data Stuff --------------------------------------------------------------

panela <- tibble(
    
    estimate = c(.831, 0.624, 0.149, 0.177, 1.339, 1.095, 0.102, 0.137),
    estimate_label = c("EA PGI", 
                        "EA PGI", "Father's EA", "Mother's EA",
                       "EA PGI", 
                       "EA PGI", "Father's EA", "Mother's EA"),
    se = c(0.026, 0.024, 0.009, 0.011, 0.041, 0.043, 0.010, 0.011),
    origin = c("Original", "Original", "Original", "Original", 
               "Corrected", "Corrected", "Corrected", "Corrected"),
    controls = c("Without Controls", 
                 "With Controls", "With Controls", "With Controls",
                 "Without Controls",
                 "With Controls", "With Controls", "With Controls")
    
    
)

# factoring stuff
panela <- panela %>%
    mutate(controls = factor(controls, levels = c("Without Controls", "With Controls")),
           estimate_label = factor(estimate_label, levels = c("Mother's EA", "Father's EA", "EA PGI")))

panela <- panela %>%
    mutate(est_lo = estimate - 1.96 * se,
           est_hi = estimate + 1.96 * se)



panelb <- tibble(
    estimate = c(0.095, 0.069, -0.047, 
                 0.055, 0.031, 0.068,
                 0.173, 0.062, -0.088,
                 0.108, 0.035, 0.103),
    estimate_label = c("EA PGI", "Family SES", "EA PGI X Family SES",
                       "EA PGI", "Family SES", "EA PGI X Family SES",
                       "EA PGI", "Family SES", "EA PGI X Family SES", 
                       "EA PGI", "Family SES", "EA PGI X Family SES"),
    se = c(0.008, 0.009, 0.009, 
           0.008, 0.010, 0.010,
           0.013, 0.009, 0.015,
           0.015, 0.010, 0.016),
    origin = c("Original", "Original", "Original", 
               "Original", "Original", "Original",
               "Corrected", "Corrected", "Corrected", 
               "Corrected", "Corrected", "Corrected"),
    edu = c("High School", "High School", "High School",
            "College", "College", "College",
            "High School", "High School", "High School",
            "College", "College", "College")
)


panelb <- panelb %>%
    mutate(est_lo = estimate - 1.96 * se,
           est_hi = estimate + 1.96 * se,
           edu = factor(edu, levels = c("High School", "College")))

# Plotting ----------------------------------------------------------------

line <- linesGrob(unit(1, "npc"), unit(c(.05, .95), "npc"),
                          gp = gpar(col = 'darkgrey', lwd = 2))


# ggforest ----------------------------------------------------------------
forestplot(
    df = panela,
    name = estimate_label,
    estimate = estimate,
    se = se,
    colour = origin,
    shape = origin
) +
    # facet_col(
    #     facets = ~controls,
    #     scales = "free_y",
    #     space = "free"
    # ) +
    facet_grid(controls~., scales = "free", space = "free", switch="y") +
    labs(color = NULL, shape = NULL,
         x = NULL, title = NULL) +
    scale_shape_manual(values = c(16, 17)) +
    theme(legend.background = element_rect(color = "black"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0,
                                           vjust = 0.85,
                                           hjust = 1)) +
    annotation_custom(line, xmin = -0.3, xmax = -0.3) +
    coord_cartesian(clip = "off")

ggsave("forest_plota.png",
       height = 6, width = 10)

forestplot(
    df = panelb,
    name = estimate_label,
    estimate = estimate,
    se = se,
    colour = origin,
    shape = origin
) +
    # facet_col(
    #     facets = ~edu,
    #     scales = "free_y",
    #     space = "free"
    # ) +
    facet_grid(edu~., scales = "free", space = "free", switch="y") +
    labs(color = NULL,  shape = NULL,
         x = NULL, title = NULL) +
    scale_shape_manual(values = c(16, 17)) +
    theme(legend.background = element_rect(color = "black"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0,
                                           vjust = 0.85,
                                           hjust = 1))+
    annotation_custom(line, xmin = -0.225, xmax = -0.225) +
    coord_cartesian(clip = "off")


ggsave("forest_plotb.png",
       height = 6, width = 10)


forestplot(
    df = panela,
    name = estimate_label,
    estimate = estimate,
    se = se,
    colour = origin,
    shape = origin
) +
    facet_col(
        facets = ~controls,
        scales = "free_y",
        space = "free"
    ) +
    labs(color = NULL, shape = NULL,
         x = NULL, title = NULL) +
    scale_shape_manual(values = c(16, 17)) +
    theme(legend.background = element_rect(color = "black"))

ggsave("forest_plota_alt.png",
       height = 6, width = 8)

forestplot(
    df = panelb,
    name = estimate_label,
    estimate = estimate,
    se = se,
    colour = origin,
    shape = origin
) +
    facet_col(
        facets = ~edu,
        scales = "free_y",
        space = "free"
    ) +
    labs(color = NULL,  shape = NULL,
         x = NULL, title = NULL) +
    scale_shape_manual(values = c(16, 17)) +
    theme(legend.background = element_rect(color = "black"))


ggsave("forest_plotb_alt.png",
       height = 6, width = 8)
