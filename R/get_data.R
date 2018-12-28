#' @title make_portal_block
#' @description Create a data.frame from the Portal rodent data with specified
#'   arguments, each row corresponds to a newmmoonnumber, and missing data are
#'   interpolated.
#' @param filter_q the numerical quantile by which to filter species. Only
#'   species that are present at least `filter_q` fraction of the time are
#'   included; default (NULL) keeps all species.
#' @inheritParams portalr::get_rodent_data
#'
#' @return a data.frame with columns for `censusdate`, and each species
#'
#' @export
make_portal_block <- function(filter_q = NULL, output = "abundance",
                              plots = c(2, 4, 8, 11, 12, 14, 17, 22), ...)
{
    # options here are:
    #   time = "all" (allow us to do correct interpolation and accouting for seasonality)
    #   level = "plot" (allow us to pull out abundances on the plots we want)
    #   effort = TRUE (so we can check effort)
    #   na_drop = TRUE (ignore periods where sampling did not occur)
    raw_rodent_data <- portalr::get_rodent_data(time = "all",
                                                plots = plots,
                                                effort = TRUE,
                                                na_drop = TRUE,
                                                output = output,
                                                ...)

    # summarize by each newmmonnumber, and for only the control plots we want
    block <- raw_rodent_data %>%
        dplyr::filter(censusdate < "2015-04-18") %>%
        dplyr::select(-period)

    # check that effort is equal across samples
    stopifnot(length(unique(block$ntraps)) == 1)

    if (!is.null(filter_q))
    {
        species_list <- block %>%
            tidyr::gather(species, abundance, BA:SO) %>%
            dplyr::group_by(species) %>%
            dplyr::summarize(quantile_q = quantile(abundance, 1 - filter_q)) %>%
            dplyr::filter(quantile_q > 0) %>%
            dplyr::pull(species)

        block <- block %>%
            dplyr::select(c("newmoonnumber", "censusdate", "ntraps", "nplots", species_list))
    }

    # add in NAs for unsampled newmoonnumbers and interpolate
    block <- block %>%
        tidyr::complete(newmoonnumber = tidyr::full_seq(newmoonnumber, 1), fill = list(NA)) %>%
        dplyr::mutate_at(dplyr::vars(-newmoonnumber, -ntraps), forecast::na.interp) %>%
        dplyr::mutate_at(dplyr::vars(-newmoonnumber, -ntraps), as.numeric) %>%
        dplyr::mutate(censusdate = as.Date(as.numeric(censusdate), origin = "1970-01-01")) %>%
        dplyr::select(-newmoonnumber, -ntraps, -nplots)

    return(block)
}

#' @title make_surrogate_annual_spline
#' @description Generate surrogate time series by computing a mean seasonal
#'   trend for each year, and shuffling the residuals. This differs from
#'   \code{\link[rEDM]{make_surrogate_seasonal}} in that the data are not
#'   sampled uniformly at the same time each year. Thus, we also ask for the
#'   `day_of_year` as input, in order to compute the mean seasonal trend.
#' @param day_of_year a vector of the numerical day of year, e.g. January 1 = 1,
#'   January 2 = 2, December 31 = 365 (leap years are a bit funny, but I don't
#'   think it should have a large effect)
#' @inheritParams rEDM::make_surrogate_data
#'
#' @return A matrix where each column is a separate surrogate with the same
#'   length as `ts`.
#'
#' @export
make_surrogate_annual_spline <- function(day_of_year, ts, num_surr = 100)
{
    # filter out NA first
    n <- length(day_of_year)
    idx <- which(is.finite(day_of_year) & is.finite(y))
    day_of_year <- day_of_year[idx]
    ts <- ts[idx]

    xx <- c(day_of_year - 365, day_of_year, day_of_year + 365)
    yy <- c(ts, ts, ts)
    seasonal_F <- smooth.spline(xx, yy)
    seasonal_cyc <- predict(seasonal_F, day_of_year)$y
    seasonal_resid <- ts - seasonal_cyc

    vals <- matrix(unlist(
        lapply(seq(num_surr), function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        })
    ), ncol = num_surr)

    out <- matrix(NA, nrow = n, ncol = num_surr)
    out[idx, ] <- vals
}
