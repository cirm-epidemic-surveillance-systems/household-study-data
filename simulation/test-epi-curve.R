
library(tidyverse)

options(scipen = 100)

# dummy normalized alpha(t)
x <- seq(0, 365, 1)
y <- dlnorm(x, meanlog = 4.5, sdlog = 0.4)
df <- data.frame(date = x, inst = y)
write.table(df, "simulation/epi-curve.csv", row.names = FALSE, quote = FALSE, col.names=F, sep = ",")

ggplot(df, aes(x = date, y = inst)) +
  geom_line()

model_choice <- "simulation/simulation-model"
hhmodel <- housepi::HouseholdModel$new(model_choice)

hhmodel$load_param_list(
  delta_t = 0.04,
  alpha = 0.1,
  beta = 0.2,
  gamma = 1,
  size_ref = 4,
  ctct_child2child = 1.5,
  ctct_child2adult = 1.2,
  ctct_child2senior = 0.5,
  ctct_adult2child = 1.2,
  ctct_adult2senior = 0.75,
  ctct_senior2child = 0.5,
  ctct_senior2adult = 0.75,
  ctct_senior2senior = 0.75
)

hhmodel$load_augm_list(
  infec_date = 1000,
  onset_date = 1000,
  recov_date = 1000,
  incub_period = list(model="reparam_gamma", args = c(mean = 3, var = 0.75)),
  ip_mean = 2,
  ip_var = 1
)

hhmodel$dataset <- housepi::survey
simulated <- hhmodel$simulate(progress = TRUE, seed = 3443)

head(simulated)
#write.csv(simulated, "dummy-output.csv", row.names = F, quote = F)
