#####################
### RUN INFERENCE ###
#####################

library(tidyverse)

bayesplot::bayesplot_theme_set(ggplot2::theme_gray())
bayesplot::color_scheme_set(scheme = "brightblue")

# n. of openMP threads over households (to use on the cluster)
n_threads <- 24

# working directory (to use on the cluster)
#work_dir <- "household-study-data"
work_dir <- "."

out_dir <- file.path(work_dir, "inference", "tmp")
if (!dir.exists(out_dir)) dir.create(out_dir)

############
### DATA ###
############

# load household survey collected according to perfect protocol
dataset <- read_csv(file.path(work_dir, "data", "simulated_perfect_survey.csv")) |>
  select(-id, -house_id) |>
  rename(id = id_num, house_id = house_id_num) |>
  # extract only dates of positive tests
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
  arrange(house_id)

# FIXME: check that is_symp & !is.na(onset_date)

# ASSUMPTIONS:
# 1. all onsets are due to current outbreak
# 2. there are no false positives or false negatives
dataset <- dataset |>
  mutate(
    # extract date of first positive information
    first_pos_date = ifelse(is_symp, onset_date, first_test_date),
    # assign a status
    status = ifelse(is_symp, 1, ifelse(is.na(first_test_date), 0, 2))
  ) |>
  select(-onset_date, -first_test_date, -is_symp)

# FIXME: no asymptomatic!
dataset |> pull(status) |> unique()

#############
### MODEL ###
#############

# call the custom inference model
model_choice <- file.path(work_dir, "inference", "inference-model")
hhmodel <- housepi::HouseholdModel$new(model_choice)

###########################
### STANDARD PARAMETERS ###
###########################

# associate parameters for inference
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

###########################
### AUGMENTED VARIABLES ###
###########################

hhmodel$load_augm_list(
  # first_pos_date of Symptomatic : onset date
  # first_pos_date of Asymptomatic : earliest positive test
  # note: no augmentation, read it from data
  first_pos_date = list(
    value = "first_pos_date@D",
    discr = TRUE
  ),
  # delay_period of Symptomatic : incubation period
  # incubation period according to "independent" study (~)
  delay_period = list(
    model = "gamma_meanvar", discr = TRUE,
    init_rand = TRUE,
    args = c(mean = 3, var = 1), 
    status = 1, rule = "draw"
  ),
  # delay_period of Asymptomatic : no a prior
  # max according to the testing protocol
  delay_period = list(
    model = "unif", discr = TRUE,
    init_rand = TRUE,
    args = c(min = 0, max = 3), 
    status = 2, rule = "draw"
  )
)

###################
### GROUP RATES ###
###################

# associate dataset
hhmodel$dataset <- dataset |>
  mutate(age_group = ifelse(age <= 12, "child", "adult"))

# associate group rates
hhmodel$load_group_rate_list(
  community = list(by = "age_group"),
  infectivity = list(by = "age_group"),
  susceptibility = list(by = "age_group")
)
# show current group rate definitions
hhmodel$group_rate_tables

##########################
### INFECTIVIY PROFILE ###
##########################

# generic infectivity profile
hhmodel$infec_profile_defs <- list(
  name = "gamma_meanvar", discr = TRUE,
  args = c(mean = "prof_mean@P", var = "prof_var@P")
)
# append infectivity profile arguments as parameters
hhmodel$load_param_list(
  prof_mean = list(prior = "unif", args = c(min = 0.1, max = 6), rule = "lnorm"),
  prof_var = list(prior = "unif", args = c(min = 0.01, max = 6), rule = "lnorm"),
  append = TRUE
)

##########################################
### INFERENCE MCMC METROPOLIS-HASTINGS ###
##########################################

# perform inference
hhmodel$infer(n_chains = 4, sample_size = 2500, n_adapt = 2000,
              n_burnin = 1000, thin = 5, n_threads = n_threads,
              seed = 3731, progress = 10)

# show mcmc traces
gtrace <- bayesplot::mcmc_trace(hhmodel$traces)
ggsave(file.path(out_dir,"mcmc-traces.pdf"), 
       plot = gtrace,
       width = 10, height = 6, dpi = 300, bg = "transparent")

# extract estimates
estimates <- hhmodel$estimates
results <- data.frame(param = row.names(estimates), 
                      estimate = estimates$mean,
                      cri_lower = estimates$cri_lower, 
                      cri_upper = estimates$cri_upper)
write.csv(results,
          file.path(out_dir, "mcmc-estimates.csv"), 
          row.names = FALSE, quote = FALSE)

################
### MEASURES ###
################

# plot infectivity profile
prof_mean <- results[results$param == "prof_mean", "estimate"]
prof_var <- results[results$param == "prof_var", "estimate"]

hhblank <- housepi::HouseholdModel$new("blank")
hhblank$add_infec_profile(name = "gamma_meanvar", discr = TRUE,
                          args = c(mean = prof_mean, var = prof_var))
infec <- hhblank$showcase_infec_profile()

ggplot(infec$profile, aes(x = gentime)) +
  geom_bar(aes(y = inst), stat = "identity", fill = "coral") +
  labs(
    x = "Generation Delay", y = "Daily", 
    title = sprintf("Infectivity Profile [mean=%.1f, sd=%.1f]",
                    infec$mean, sqrt(infec$var))
  ) +
  coord_cartesian(xlim = c(NA, 15))

# save infectivity profile
ggsave(file.path(out_dir, "infec_prof.pdf"),
       width = 6, height = 4, dpi = 300, bg = "transparent")
