
# Read DHS data
read_dhs_data <- function(path) {
  haven::read_dta(path) |>
    dplyr::transmute(
      house_id = hhid,
      id = NA,
      size = hv009,
      age = hv105) |>
    dplyr::group_by(house_id) |>
    dplyr::mutate(
      house_id = dplyr::cur_group_id(),
      id = dplyr::row_number())
}


# Bootstrap n household structures
sample_households <- function(n, data, max_size = Inf) {

  list_hh <- dplyr::filter(data, size <= max_size) |>
    group_by(house_id) |>
    group_split()

  idx_hh <- sample.int(length(list_hh), n, replace = TRUE)

  ret <- list_hh[idx_hh]

  ret |>
    bind_rows(.id = "new_hh_id") |>
    dplyr::transmute(id = id,
                     house_id = new_hh_id,
                     age = age,
                     size = size)
}

# Plot household structure data set
plot_households <- function(data) {

  n_hh <- max(as.integer(data$house_id))
  n <- nrow(data)

  g_age <- data |>
    ggplot(aes(x = age)) +
    geom_histogram() +
    labs(x = "Age", y = "Frequency") +
    theme_bw()

  g_size <- data |>
    ggplot(aes(x = size)) +
    geom_bar() +
    labs(x = "Household size", y = "Frequency") +
    theme_bw()


  g_age_by_size <- data |>
    ggplot(aes(x = age, group = size)) +
    geom_density_ridges(aes(y = size)) +
    labs(y = "Household size", x = "Age") +
    theme_bw()

  g <- ((g_size / g_age) | g_age_by_size) +
    plot_annotation(tag_levels = "A",
                    title = sprintf("Data contains %s individuals within %s households",
                                    prettyNum(n, big.mark = ","),
                                    prettyNum(n_hh, big.mark = ",")))
  g
}
