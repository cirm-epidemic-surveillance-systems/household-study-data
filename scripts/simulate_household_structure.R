
library(rdhs)
library(dplyr)
library(haven)
library(ggplot2)
library(ggridges)
library(patchwork)

source("R/household_structure.R")

# Read in Model DHS data
dhs_data <- read_dhs_data("data/dhs_model_data.dta")

# Take a bootstrapped sample
data_sample <- sample_households(100, dhs_data)

# Plot full DHS data set
g_dhs <- plot_households(dhs_data)
ggsave("figures/household_structure_dhs.png", g_dhs)


# Plot sampled data set
g_sample <- plot_households(data_sample)
ggsave("figures/household_structure_sample.png", g_sample)
