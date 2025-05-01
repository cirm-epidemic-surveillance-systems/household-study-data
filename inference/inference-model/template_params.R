list(
  alpha = list(
    value = 0.0001, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", scale = 1
  ),
  beta = list(
    value = 0.07, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", scale = 1
  ),
  gamma = list(
    value = 0.9, prior = "unif", args = list(min = 0, max = 4),
    rule = "norm", scale = 1
  ),
  size_ref = list(value = 4),
  study_period = list(value = 45)
)
