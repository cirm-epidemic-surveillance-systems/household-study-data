
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
plot_households(dhs_data)

# Plot sampled data set
plot_households(data_sample)
