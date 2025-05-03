
library(tidyverse)

n_threads <- 36

path <- "household-study-data"
#path <- "."

bayesplot::bayesplot_theme_set(ggplot2::theme_gray())
bayesplot::color_scheme_set(scheme = "brightblue")

## dummy dataset to play around
# proba_asym <- 0.1
# dataset <- read_csv("data/dummy-output.csv") |>
#   mutate(
#     onset_date = floor(onset_date + (onset_date %% 1)),
#     status = if_else(status == 1 & runif(n()) < !!proba_asym, 2L, status),
#     first_pos_date = if_else(status == 2, infec_date + runif(1, min = 2, max = 10), onset_date)
#   ) |>
#   select(-onset_date, -incub_period, -infec_date)
# dataset

# FIXME: missing age!
dataset <- read_csv(file.path(path,"data/simulated_perfect_survey.csv")) |>
  select(-id, -house_id) |>
  rename(id = id_num, house_id = house_id_num) |>
  # date of positive tests
  mutate(pos_test_date = ifelse(isTRUE(test_outcome), t, NA)) |>
  select(-t) |>
  group_by(id) |>
  summarise(
    house_id = unique(house_id),
    age = unique(age),
    onset_date = unique(onset_date),
    is_symp = unique(is_symp),
    # first positive test
    first_pos_test_date = ifelse(all(is.na(pos_test_date)), 
                                 NA, min(pos_test_date, na.rm = TRUE))
  ) |>
  arrange(house_id) |>
  mutate(
    # first positive information
    first_pos_date = ifelse(
      is.na(first_pos_test_date) & is_symp,
      onset_date,
      ifelse(is.na(onset_date), 
             first_pos_test_date, min(onset_date, first_pos_test_date))
    ),
    status = as.integer(!is.na(first_pos_date))
  ) |>
  select(-onset_date, -first_pos_test_date, -is_symp)

# FIXME: no asymptomatic!
dataset |> pull(status) |> unique()

model_choice <- file.path(path, "inference/inference-model")
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
  study_period = 365 # FIXME: study period according to protocol
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
    # args according to "independent" study (~)
    args = c(mean = 3, var = 1), 
    status = 1, rule = "draw"
  ),
  # delay_period of Asymptomatic : no a prior
  delay_period = list(
    model = "unif", discrete = TRUE,
    init_rand = TRUE,
    # max according to the testing protocol
    args = c(min = 0, max = 3), 
    status = 2, rule = "draw"
  )
)

dataset <- dataset |>
  mutate(age_group = ifelse(age <= 12, "child", "adult"))

hhmodel$dataset <- dataset

# associate group rates
hhmodel$load_group_rate_list(
  community = list(by = "age_group"),
  infectivity = list(by = "age_group"),
  susceptibility = list(by = "age_group")
)
hhmodel$group_rate_tables

# infectivity profile
hhmodel$infec_profile_defs <- list(name = "reparam_gamma",
                                   discrete = TRUE,
                                   args = c(mean = "prof_mean", var = "prof_var"))
hhmodel$load_param_list(
  prof_mean = list(prior = "unif", args = c(min = 0.1, max = 6), rule = "lnorm"),
  prof_var = list(prior = "unif", args = c(min = 0.01, max = 6), rule = "lnorm"),
  append = TRUE
)

# perform inference
hhmodel$infer(n_chains = 4, n_adapt = 2000, 
              sample_size = 2500,
              n_burnin = 1000, thin = 5, n_threads = n_threads,
              seed = 3731, progress = 10)

# show mcmc traces
gtrace <- bayesplot::mcmc_trace(hhmodel$traces)
ggsave(file.path(path,"inference/mcmc-traces.pdf"), 
      plot = gtrace,
      width = 10, height = 6, dpi = 300, bg = "transparent")

# extract results
x <- hhmodel$traces_results
results <- data.frame(param = row.names(x), estimate = x$mean,
                      cri_lower = x$cri_lower, cri_upper = x$cri_upper)
write.csv(results,
          file.path(path, "inference/mcmc-estimates.csv"), 
          row.names = FALSE, quote = FALSE)