
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
sample_households <- function(n, data) {
  n_hh <- max(data$house_id)
  idx_hh <- sample.int(n_hh, n)

  data |>
    dplyr::filter(house_id %in% idx_hh) |>
    dplyr::group_by(house_id) |>
    dplyr::mutate(house_id = dplyr::cur_group_id())
}

# Plot household structure data set
plot_households <- function(data) {

  n_hh <- max(data$house_id)
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
