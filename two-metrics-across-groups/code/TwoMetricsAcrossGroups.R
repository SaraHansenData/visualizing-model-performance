# Visualizing Model Performance: Two metrics across groups
## Sara Hansen, hanse2s
## Modified 2024-04-18

library(tidyverse) # data processing, visualization
library(ggthemes) # all-in-one themes
library(ggrepel) # neat labeling
library(plotrix) # standard error calculation
library(grid) # annotation text customization
library(ggimage) # custom shapes for scatterplot
library(magick) # read image for custom annotation

options(max.print = 1e9)
options(dplyr.print_max = 1e9)

# Note you will need to use ggsave() on plots to see their actual appearance,
# the built-in plot viewer (if you are using RStudio) is not accurate

# # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# In this fictional example, we built predictive models for Sulu pygmy woodpecker
# Now, we want to plot the model performance metrics to compare them

# In this script, we'll visualize two model performance metrics across three groups
# We'll create a few variants of nice plots to suit different questions

# We need three distinct colors
# https://palett.es/ provides a great randomizer with accessibility checker

# A few colors have been picked out already, but change them here if you like
col1 <- "#a12441"
col2 <- "#c05f2e"
col3 <- "#263d45"

# # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Import data
dat <- read.delim("data/model-performance-data.txt", header = TRUE,
                  sep = "\t", quote = "", fileEncoding = "UTF-8")

# Explore data
glimpse(dat)
head(dat)
summary(dat)

# For each model, we have:
  # ModelID (unique identifier)
  # Site (Lowland forest, Mangrove forest, or Cropland)
  # Covariates (Climate, Human, and Biotic and all possible combinations)
  # Precision (proportion of positive predictions that are actually positive)
  # Recall (proportion of actual positives that were predicted as positive)
  # F1 Score (metric that incorporates and balances Precision and Recall)

# The question we want to answer is:
  # How do Precision and Recall compare based on the site and the covariates used in the model?

# Plot 1: A simple scatterplot, colored by group
dat %>%
  ggplot(aes(x = Precision, y = Recall, color = Site)) +
  geom_point(size = 2) +
  xlim(0.46, 0.99) + # customize axis limits for cleaner look
  ylim(0.46, 0.99) + # customize axis limits for cleaner look
  scale_color_manual(values = c(col1, col2, col3)) + # use our own colors
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"), # add gridlines
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))
# This shows us all the models' Precision and Recall, but we can't tell a lot about the groups
#ggsave("Plot1.pdf", width = 7, height = 7, units = "in")

# Plot 2: A labeled scatterplot of model means, colored by group
dat %>%
  group_by(Site, Covariates) %>%
  summarize(meanPrecision = mean(Precision),
            meanRecall = mean(Recall)) %>% # calculate means for each model group
  ungroup() %>%
  ggplot(aes(x = meanPrecision, y = meanRecall, color = Site, label = Covariates)) +
  geom_point(size = 4) +
  ggrepel::geom_label_repel(size = 10/.pt,
                            box.padding = 0.6) + # label which point is which model group
  xlim(0.46, 0.99) +
  ylim(0.46, 0.99) +
  scale_color_manual(values = c(col1, col2, col3)) +
  labs(x = "Precision", y = "Recall") +
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))
# This shows us a mean of metrics for each covariate group, but we can't see the spread across all models
#ggsave("Plot2.pdf", width = 7, height = 7, units = "in")

# Plot 3: A labeled scatterplot of model means with error bars
dat %>%
  group_by(Site, Covariates) %>%
  summarize(meanPrecision = mean(Precision),
            meanRecall = mean(Recall),
            sePrecision = std.error(Precision),
            seRecall = std.error(Recall)) %>% # include standard error calculation
  mutate(upperPrecision = meanPrecision + sePrecision,
         lowerPrecision = meanPrecision - sePrecision,
         upperRecall = meanRecall + seRecall,
         lowerRecall = meanRecall - seRecall) %>% # upper and lower limit based on SE
  ungroup() %>%
  ggplot(aes(x = meanPrecision, y = meanRecall, color = Site, label = Covariates)) +
  geom_pointrange(aes(xmin = lowerPrecision, xmax = upperPrecision),
                  linewidth = 2/.pt, linetype = 1) + # Precision ranges
  geom_pointrange(aes(ymin = lowerRecall, ymax = upperRecall),
                  linewidth = 2/.pt, linetype = 1) + # Recall ranges
  # geom_point() is not needed when we use geom_pointran
  ggrepel::geom_label_repel(size = 10/.pt,
                            box.padding = 0.5) +
  xlim(0.55, 0.95) +
  ylim(0.55, 0.95) +
  scale_color_manual(values = c(col1, col2, col3)) +
  labs(x = "Precision", y = "Recall",
       color = "Site") +
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))
# We can now see means and standard error, but it looks messy with all the covariate labels
#ggsave("Plot3.pdf", width = 7, height = 7, units = "in")

# Plot 4: An unlabeled scatterplot of model means with error bars
# This is the same as Plot 3 but removing labeling and using different shapes instead
dat %>%
  group_by(Site, Covariates) %>%
  summarize(meanPrecision = mean(Precision),
            meanRecall = mean(Recall),
            sePrecision = std.error(Precision),
            seRecall = std.error(Recall)) %>% # include standard error calculation
  mutate(upperPrecision = meanPrecision + sePrecision,
         lowerPrecision = meanPrecision - sePrecision,
         upperRecall = meanRecall + seRecall,
         lowerRecall = meanRecall - seRecall) %>% # upper and lower limit based on SE
  ungroup() %>%
  ggplot(aes(x = meanPrecision, y = meanRecall, color = Site, shape = Covariates)) +
  geom_pointrange(aes(xmin = lowerPrecision, xmax = upperPrecision),
                  linewidth = 2/.pt, linetype = 1, size = 0.65) + # Precision ranges
  geom_pointrange(aes(ymin = lowerRecall, ymax = upperRecall),
                  linewidth = 2/.pt, linetype = 1, size = 0.65) + # Recall ranges
  xlim(0.55, 0.95) +
  ylim(0.55, 0.95) +
  scale_color_manual(values = c(col1, col2, col3)) +
  labs(x = "Precision", y = "Recall",
       color = "Site") +
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.box = "vertical") # stack the two legends
# This is better, but the covariate shapes aren't intuitive
#ggsave("Plot4.pdf", width = 7.5, height = 7, units = "in")

# Plot 5: An unlabeled scatterplot of model means with error bars
# and iconic shapes for the covariates
# This is the same as Plot 4 with custom shapes, which we'll download from GitHub
dat %>%
  group_by(Site, Covariates) %>%
  summarize(meanPrecision = mean(Precision),
            meanRecall = mean(Recall),
            sePrecision = std.error(Precision),
            seRecall = std.error(Recall)) %>%
  mutate(upperPrecision = meanPrecision + sePrecision,
         lowerPrecision = meanPrecision - sePrecision,
         upperRecall = meanRecall + seRecall,
         lowerRecall = meanRecall - seRecall) %>%
  mutate(image = case_when(
    Covariates == "Biotic" ~  "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Biotic.png?raw=true",
    Covariates == "Biotic + Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Biotic%2BHuman.png?raw=true",
    Covariates == "Climate" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate.png?raw=true",
    Covariates == "Climate + Biotic" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate%2BBiotic.png?raw=true",
    Covariates == "Climate + Biotic + Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate%2BBiotic%2BHuman.png?raw=true",
    Covariates == "Climate + Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate%2BHuman.png?raw=true",
    Covariates == "Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Human.png?raw=true")) %>%
  ungroup() %>%
  ggplot(aes(x = meanPrecision, y = meanRecall)) +
  geom_point(aes(color = Site), size = 8, shape = 16) +
  # geom_point() here creates a color-coded outline around the images
  geom_pointrange(aes(xmin = lowerPrecision, xmax = upperPrecision,
                      color = Site),
                  linewidth = 2/.pt, linetype = 1) + # Precision ranges
  geom_pointrange(aes(ymin = lowerRecall, ymax = upperRecall,
                      color = Site),
                  linewidth = 2/.pt, linetype = 1) + # Recall ranges
  geom_image(aes(image = image), size = 0.035) +
  coord_cartesian(xlim = c(0.55, 0.95),
                  ylim = c(0.55, 0.95),
                  clip = "off") +
  # we don't need xlim and ylim because we are using coord_cartesian so that we can annotate outside of the plot
  scale_color_manual(values = c(col1, col2, col3)) +
  labs(x = "Precision", y = "Recall",
       color = "Site") +
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
        plot.margin = margin(1, 1, 2, 1, unit = "cm"), # making space for the annotation
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(color = NA)) + # removing Site legend border
  # now we'll build our own legend for the covariate icons
  annotation_custom(grid::textGrob(expression(bold("Covariates")), gp = gpar(fontsize = 12)), 0.64, 0.64, 0.405, 0.405) +
  annotation_custom(grid::rasterGrob(image_read("https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Biotic.png?raw=true"),
                                     interpolate = TRUE),
                    0.68, 0.71, 0.39, 0.42) +
  annotation_custom(grid::textGrob("Biotic", gp = gpar(fontsize = 10)), 0.72, 0.72, 0.405, 0.405) +
  annotation_custom(grid::rasterGrob(image_read("https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate.png?raw=true"),
                                     interpolate = TRUE),
                    0.745, 0.775, 0.39, 0.42) +
  annotation_custom(grid::textGrob("Climate", gp = gpar(fontsize = 10)), 0.79, 0.79, 0.405, 0.405) +
  annotation_custom(grid::rasterGrob(image_read("https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Human.png?raw=true"),
                                     interpolate = TRUE),
                    0.82, 0.85, 0.39, 0.42) +
  annotation_custom(grid::textGrob("Human", gp = gpar(fontsize = 10)), 0.865, 0.865, 0.405, 0.405)
# Note we placed some of the aes() arguments in different places here
# to avoid the issue of coloring over the images
# We also annotated over the plot to get the icons onto the legend 
# because ggimage does not yet allow for custom legend keys
#ggsave("Plot5.pdf", width = 7.5, height = 7, units = "in")


# Now what if we want to see the full range of values for a given Site's models?
# We'll add convex hulls to see the distribution of Precision and Recall overall
dat2 <- dat %>%
  select(Site, x = Precision, y = Recall) %>%
  group_by(Site)

find_hull <- function(points) points[chull(points$x, points$y), ]
hulls <- plyr::ddply(dat2, "Site", find_hull)

# We can add these convex hulls to any of the above plots
# For example
# Plot 6: Plot 1 plus convex hulls
dat %>%
  ggplot(aes(x = Precision, y = Recall, color = Site)) +
  geom_point(size = 2) +
  xlim(0.46, 0.99) + # customize axis limits for cleaner look
  ylim(0.46, 0.99) + # customize axis limits for cleaner look
  scale_color_manual(values = c(col1, col2, col3)) + # use our own colors
  geom_polygon(data = hulls, aes(x = x, y = y, fill = Site, label = NULL),
               alpha = 0.1, linetype = 0) +
  scale_fill_manual(values = c(col1, col2, col3)) +
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"), # add gridlines
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))
# We don't gain a lot of information, but we get a sense of the overlap between Site models
#ggsave("Plot6.pdf", width = 7, height = 7, units = "in")

# Plot 7: Bringing it all together into one informative plot
dat %>%
  group_by(Site, Covariates) %>%
  summarize(meanPrecision = mean(Precision),
            meanRecall = mean(Recall),
            sePrecision = std.error(Precision),
            seRecall = std.error(Recall)) %>%
  mutate(upperPrecision = meanPrecision + sePrecision,
         lowerPrecision = meanPrecision - sePrecision,
         upperRecall = meanRecall + seRecall,
         lowerRecall = meanRecall - seRecall) %>%
  mutate(image = case_when(
    Covariates == "Biotic" ~  "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Biotic.png?raw=true",
    Covariates == "Biotic + Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Biotic%2BHuman.png?raw=true",
    Covariates == "Climate" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate.png?raw=true",
    Covariates == "Climate + Biotic" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate%2BBiotic.png?raw=true",
    Covariates == "Climate + Biotic + Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate%2BBiotic%2BHuman.png?raw=true",
    Covariates == "Climate + Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate%2BHuman.png?raw=true",
    Covariates == "Human" ~ "https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Human.png?raw=true")) %>%
  ungroup() %>%
  ggplot(aes(x = meanPrecision, y = meanRecall)) +
  geom_point(aes(color = Site), size = 8, shape = 16) +
  geom_pointrange(aes(xmin = lowerPrecision, xmax = upperPrecision,
                      color = Site),
                  linewidth = 2/.pt, linetype = 1) +
  geom_pointrange(aes(ymin = lowerRecall, ymax = upperRecall,
                      color = Site),
                  linewidth = 2/.pt, linetype = 1) +
  geom_image(aes(image = image), size = 0.035) +
  coord_cartesian(xlim = c(0.51, 0.94),
                  ylim = c(0.46, 0.99),
                  clip = "off") +
  scale_color_manual(values = c(col1, col2, col3)) +
  labs(x = "Precision", y = "Recall",
       color = "Site") +
  geom_polygon(data = hulls, aes(x = x, y = y, fill = Site, label = NULL),
               alpha = 0.1, linetype = 0) +
  scale_fill_manual(values = c(col1, col2, col3)) +
  theme_clean() +
  theme(panel.grid.major.x = element_line(colour = "gray", linetype = "dotted"),
        plot.margin = margin(1, 1, 1.8, 1, unit = "cm"),
        legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(color = NA)) + 
  annotation_custom(grid::textGrob(expression(bold("Covariates")), gp = gpar(fontsize = 12)), 0.615, 0.615, 0.28, 0.28) +
  annotation_custom(grid::rasterGrob(image_read("https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Biotic.png?raw=true"),
                                     interpolate = TRUE),
                    0.655, 0.685, 0.265, 0.295) +
  annotation_custom(grid::textGrob("Biotic", gp = gpar(fontsize = 10)), 0.695, 0.695, 0.28, 0.28) +
  annotation_custom(grid::rasterGrob(image_read("https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Climate.png?raw=true"),
                                     interpolate = TRUE),
                    0.72, 0.75, 0.265, 0.295) +
  annotation_custom(grid::textGrob("Climate", gp = gpar(fontsize = 10)), 0.765, 0.765, 0.28, 0.28) +
  annotation_custom(grid::rasterGrob(image_read("https://github.com/hanse2s/visualizing-model-performance/blob/main/custom-shapes/Human.png?raw=true"),
                                     interpolate = TRUE),
                    0.795, 0.825, 0.265, 0.295) +
  annotation_custom(grid::textGrob("Human", gp = gpar(fontsize = 10)), 0.84, 0.84, 0.28, 0.28)
#ggsave("Plot7.pdf", width = 7.5, height = 7, units = "in")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

rm(list = ls())



