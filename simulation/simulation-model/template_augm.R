list(
  onset_date = list(
    value = "onset_date@D",
    rule = "norm", scale = 0.5,
    from = "onset_date@D"
  ),
  incub_period = list(
    init_rand = TRUE, model = "gamma",
    args = c(shape = 2, scale = 1),
    max = 14, rule = "gibbs"
  )
)
