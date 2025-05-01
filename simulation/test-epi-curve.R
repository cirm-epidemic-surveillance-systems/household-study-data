
library(tidyverse)
library(patchwork)
library(ggridges)

source("R/household_structure.R")

options(scipen = 100)

# dummy normalized alpha(t)
x <- seq(0, 365, 1)
y <- dlnorm(x, meanlog = 4.5, sdlog = 0.4)
df <- data.frame(date = x, inst = y)
write.table(df, "simulation/epi-curve.csv", row.names = FALSE, quote = FALSE,
            col.names = F, sep = ",")

ggplot(df, aes(x = date, y = inst)) +
  geom_line()

model_choice <- "simulation/simulation-model"
hhmodel <- housepi::HouseholdModel$new(model_choice)

hhmodel$load_param_list(
  delta_t = 0.04,
  alpha = 0.05,
  beta = 0.4,
  gamma = 1,
  size_ref = 4,
  ctct_child2child = 1.5,
  ctct_child2adult = 1.2,
  ctct_child2senior = 0.5,
  ctct_adult2child = 1.2,
  ctct_adult2senior = 0.75,
  ctct_senior2child = 0.5,
  ctct_senior2adult = 0.75,
  ctct_senior2senior = 0.75,
  com_child = 1,
  com_senior = 0.5,
  inf_child = 1,
  inf_senior = 1,
  sus_child = 1,
  sus_senior = 0.4
)

# In gamma dist mean = shape / rate
# var = mean / rate -> rate = mean / var
# shape = rate * mean = mean^2 / var
# so if shape > 1 we need mean^2 > var

hhmodel$load_augm_list(
  infec_date = 1000,
  onset_date = 1000,
  recov_date = 1000,
  incub_period = list(model="reparam_gamma", args = c(mean = 3, var = 0.75)),
  ip_mean = list(model = "unif", args = c(min = 2, max = 4)),
  ip_var = list(model = "unif", args = c(min = 0.5, max = 2))
)

# bootstrap dhs household data
dhs_data <- read_dhs_data("data/dhs_model_data.dta")
n_hh <- 1e3
hhmodel$dataset <- sample_households(n_hh, dhs_data)

plot_households(hhmodel$dataset)

# Run infection model
simulated <- hhmodel$simulate(progress = TRUE, seed = 3443)

write.csv(simulated, "simulation/outputs/simulated_infections.csv",
          row.names = F, quote = F)

# Exploratory plots of infection process
infected <- simulated |>
  filter(status == 1)

g <- infected |>
  group_by(floor(infec_date)) |>
  count() |>
  ungroup() |>
  mutate(prop = n / sum(n)) |>
  ggplot(aes(x = `floor(infec_date)`, y = prop)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  geom_line(data = df, aes(x = date, y = inst)) +
  labs(x = "Day", y = "Proportion of infections occuring on each day")
ggsave("simulation/outputs/fig_daily_infections_vs_external_foi.png", g)

final_size_individual <- simulated |>
  summarise(final_size = sum(status),
            total_pop = n(),
            seroprevalence = final_size / total_pop)

final_size_hh <- simulated |>
  group_by(house_id) |>
  mutate(hh_inf = sum(status)) |>
  ungroup() |>
  summarise(n = sum(hh_inf > 0),
            N = n(),
            prop = n / N)

g <- infected |>
  mutate(age_group = cut(age, c(0, 16, 65, 120), right = FALSE)) |>
  group_by(floor(infec_date), age_group) |>
  count() |>
  ggplot(aes(x = `floor(infec_date)`, y = n, fill = age_group)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  theme_bw() +
  labs(x = "Day", y = "Daily number of infections",
       title = "")
ggsave("simulation/outputs/fig_daily_infections_by_age_group.png", g)


who_stats <- simulated |>
  group_by(house_id) |>
  summarise(size = n(), n = sum(status)) |>
  filter(n > 0) |> # constrain to households with at least one case
  mutate(SIR = (n - 1) / (size - 1))

who_stats |>
  ggplot(aes(x = SIR)) +
  geom_histogram() +
  theme_bw()

g_a <- who_stats |>
  group_by(size) |>
  summarise(mean = mean(SIR)) |>
  ggplot(aes(x = size, y = mean)) +
  geom_point() +
  theme_bw() +
  labs(x = "Household size", y = "Mean SIR")

g_b <- who_stats |>
  ggplot(aes(x = size)) +
  geom_bar() +
  theme_bw() +
  labs(x = "Household size", y = "Frequency")

g <- g_a / g_b +
  plot_annotation(tag_levels = "A")
ggsave("simulation/outputs/fig_mean_SIR.png", g)
