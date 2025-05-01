library(tidyverse)

options(scipen = 100)

model_choice <- "inference/inference-model"
hhmodel <- housepi::HouseholdModel$new(model_choice)

hhmodel$load_param_list(
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
  # first_pos_date of Symptomatic : onset date
  # first_pos_date of Asymptomatic : earliest positive test
  # note: no augmentation, read it from data
  first_pos_date = "first_pos_date@D",
  # delay_period of Symptomatic : incubation period
  delay_period = list(
    model = "reparam_gamma",
    init_rand = TRUE,
    args = c(mean = 3, var = 0.75), # FIXME: args according to independent study
    status = 1
  ),
  # delay_period of Asymptomatic : no a prior
  delay_period = list(
    model = "unif",
    init_rand = TRUE,
    args = c(min = 0, max = 3), # FIXME: max according to protocol
    status = 2
  )
)

simulated <- hhmodel$infer()
