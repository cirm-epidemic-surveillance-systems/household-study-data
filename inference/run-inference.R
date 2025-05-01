
library(tidyverse)

# dummy dataset to play around
proba_asym <- 0.1
dataset <- read_csv("data/dummy-output.csv") |>
  mutate(
    onset_date = floor(onset_date + (onset_date %% 1)),
    status = if_else(status == 1 & runif(n()) < !!proba_asym, 2L, status),
    first_pos_date = if_else(status == 2, infec_date + runif(1, min = 2, max = 10), onset_date)
  ) |>
  select(-onset_date, -incub_period, -infec_date)
dataset

model_choice <- "inference/inference-model"
hhmodel <- housepi::HouseholdModel$new(model_choice)

hhmodel$load_param_list(
  alpha = list(
    prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", scale = 1
  ),
  beta = list(
    prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", scale = 1
  ),
  gamma = list(
    prior = "unif", args = list(min = 0, max = 4),
    rule = "norm", scale = 1
  ),
  size_ref = 4,
  study_period = 45 # FIXME: study period according to protocol
)

hhmodel$load_augm_list(
  # first_pos_date of Symptomatic : onset date
  # first_pos_date of Asymptomatic : earliest positive test
  # note: no augmentation, read it from data
  first_pos_date = list(
    value = "first_pos_date@D",
    discrete = TRUE
  ),
  # delay_period of Symptomatic : incubation period
  delay_period = list(
    model = "reparam_gamma", discrete = TRUE,
    init_rand = TRUE,
    args = c(mean = 3, var = 0.75), # FIXME: args according to independent study
    status = 1, rule = "draw"
  ),
  # delay_period of Asymptomatic : no a prior
  delay_period = list(
    model = "unif", discrete = TRUE,
    init_rand = TRUE,
    args = c(min = 0, max = 3), # FIXME: max according to protocol
    status = 2, rule = "draw"
  )
)

hhmodel$dataset <- dataset |>
  mutate(age_group = ifelse(age <= 12, "child", "adult"))

# associate group rates
hhmodel$load_group_rate_list(
  community = list(by = "age_group"),
  infectivity = list(by = "age_group"),
  susceptibility = list(by = "age_group")
)
hhmodel$group_rate_tables

# perform inference
hhmodel$infer(n_chains = 3, n_adapt = 2000, n_steps = 10000,
              n_burnin = 1000, thin = 10, seed = 3731, progress = 30)
# show mcmc traces
bayesplot::mcmc_trace(hhmodel$traces)

# extract results
hhmodel$traces_results