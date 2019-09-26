#' @title Get a data.frame of Portal rodent abundances
#' @description Create a data.frame from the Portal rodent data with specified
#'   arguments, each row corresponds to a newmmoonnumber, and missing data are
#'   interpolated.
#' @param filter_q the numerical quantile by which to filter species. Only
#'   species that are present at least `filter_q` fraction of the time are
#'   included; default (NULL) keeps all species.
#' @inheritParams portalr::summarize_rodent_data
#'
#' @return a data.frame with columns for `censusdate`, and each species
#'
#' @export
make_portal_block <- function(path = portalr::get_default_data_path(),
                              filter_q = NULL, output = "abundance",
                              plots = c(2, 4, 8, 11, 12, 14, 17, 22), ...) {
  # options here are:
  #   time = "all" (allow us to do correct interpolation and accouting for seasonality)
  #   level = "plot" (allow us to pull out abundances on the plots we want)
  #   effort = TRUE (so we can check effort)
  #   na_drop = TRUE (ignore periods where sampling did not occur)
  raw_rodent_data <- portalr::summarize_rodent_data(
    path = path,
    time = "all",
    plots = plots,
    effort = TRUE,
    na_drop = TRUE,
    output = output,
    ...
  )

  # summarize by each newmmonnumber, and for only the control plots we want
  block <- raw_rodent_data %>%
    dplyr::filter(.data$censusdate < "2015-04-18") %>%
    dplyr::select(-.data$period)

  # check that effort is equal across samples
  stopifnot(length(unique(block$ntraps)) == 1)

  if (!is.null(filter_q)) {
    species_list <- block %>%
      tidyr::gather("species", "abundance", dplyr::matches("^[A-Z]{2}$")) %>%
      dplyr::group_by(.data$species) %>%
      dplyr::summarize(quantile_q = stats::quantile(.data$abundance, 1 - filter_q)) %>%
      dplyr::filter(.data$quantile_q > 0) %>%
      dplyr::pull(.data$species)

    block <- block %>%
      dplyr::select(c("newmoonnumber", "censusdate", "ntraps", "nplots", species_list))
  }

  # add in NAs for unsampled newmoonnumbers and interpolate
  block <- block %>%
    tidyr::complete(newmoonnumber = tidyr::full_seq(.data$newmoonnumber, 1), fill = list(NA)) %>%
    dplyr::mutate_at(dplyr::vars(-.data$newmoonnumber, -.data$ntraps), forecast::na.interp) %>%
    dplyr::mutate_at(dplyr::vars(-.data$newmoonnumber, -.data$ntraps), as.numeric) %>%
    dplyr::mutate(censusdate = as.Date(as.numeric(.data$censusdate), origin = "1970-01-01")) %>%
    dplyr::select(-.data$newmoonnumber, -.data$ntraps, -.data$nplots)

  return(block)
}
