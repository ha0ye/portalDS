#' @title Run univariate analysis on each time series
#' @description Run simplex projection models for each time series in the
#'   `block` and save the output, to determine the best embedding dimension for
#'   each column.
#' @param block a data.frame containing time series for the community. Each
#'   column is a time series of abundances.
#' @param E_list the embedding dimension or range of embedding dimensions to
#'   search over.
#' @param surrogate_method which surrogate method to use:
#'   options are "annual_spline" or methods available in
#'   [rEDM::make_surrogate_data()]
#' @param num_surr number of surrogates to compute
#' @param surr_params a list of named optional arguments to be passed into the
#'   surrogate data function
#' @return A tibble with columns for the species name (taken from the original
#'   column names), the abundance time series for each species, the output from
#'   [rEDM::simplex()], the best embedding dimension, as determined by 
#'   the E that minimizes MAE, and surrogate time series
#'
#' @export
compute_simplex <- function(block, E_list = 1:10,
                            surrogate_method = "annual_spline",
                            num_surr = 100,
                            surr_params = NULL) {
  simplex_results <- block %>%
    dplyr::select(-.data$censusdate) %>%
    tidyr::gather("species", "abundance") %>%
    dplyr::group_by(.data$species) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      simplex_out =
        purrr::map(.data$data, ~ rEDM::simplex(.$abundance, E = E_list, silent = TRUE))
    ) %>%
    dplyr::mutate(best_E = purrr::map_int(
      .data$simplex_out,
      ~ dplyr::filter(., .data$mae == min(.data$mae)) %>%
        dplyr::pull(.data$E) %>%
        as.integer()
    )) %>%
    dplyr::ungroup()

  surrogate_method <- tolower(surrogate_method)
  if (surrogate_method == "twin") {
    simplex_results$surrogate_data <-
      purrr::pmap(
        dplyr::select(simplex_results, c("data", "best_E")),
        ~ do.call(
          rEDM::make_surrogate_twin,
          c(
            list(ts = ..1, dim = ..2, num_surr = num_surr),
            surr_params
          )
        )
      )
  } else if (surrogate_method == "annual_spline") {
    day_of_year <- lubridate::yday(block$censusdate)
    simplex_results$surrogate_data <-
      purrr::pmap(
        dplyr::select(simplex_results, "data"),
        ~ do.call(
          make_surrogate_annual_spline,
          c(
            list(ts = ..1, num_surr = num_surr, day_of_year = day_of_year), 
            surr_params
          )
        )
      )
  } else {
    simplex_results$surrogate_data <-
      purrr::pmap(
        dplyr::select(simplex_results, "data"),
        ~ do.call(
          rEDM::make_surrogate_data,
          c(
            list(ts = ..1, num_surr = num_surr, method = surrogate_method),
            surr_params
          )
        )
      )
  }
  return(simplex_results)
}
