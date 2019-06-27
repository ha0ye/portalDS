#' @title Run CCM on each pair of time series
#' @description `compute_ccm` runs pairwise CCM based on the simplex_output - 
#'   using the best embedding dimension from the simplex results, and computed 
#'   for both the real data and the surrogate data. The calculations run using 
#'   \code{\link[furrr]{future_pmap}}. Thus, parallelization should be set by 
#'   the user, if desired, using \code{\link[future]{plan}}, prior to running.
#' @param simplex_results the output of \code{\link{compute_simplex}}
#' @inheritParams rEDM::ccm
#' @return `compute_ccm` returns a tibble with columns for the variables that we 
#'   use in CCM, the data type (whether it's the "actual" time series or 
#'   "surrogate"), the library size, and the results from CCM
#'
#' @export
compute_ccm <- function(simplex_results,
                        lib_sizes = seq(10, 100, by = 10),
                        random_libs = TRUE, num_samples = 100,
                        replace = TRUE, RNGseed = 42,
                        silent = TRUE)
{
    ccm_func <- function(df, E)
    {
        df %>%
            rEDM::ccm(E = E, lib_sizes = lib_sizes,
                      random_libs = random_libs, num_samples = num_samples,
                      replace = replace,
                      lib_column = 1, target_column = 2,
                      RNGseed = RNGseed, silent = silent) %>%
            rEDM::ccm_means(na.rm = TRUE) %>%
            dplyr::select(lib_size, num_pred, rho, mae, rmse)
    }
    
    params <- expand.grid(from_idx = seq(NROW(simplex_results)),
                          to_idx = seq(NROW(simplex_results))) %>%
        dplyr::mutate(E = simplex_results$best_E[from_idx], 
                      from_var = simplex_results$species[from_idx], 
                      to_var = simplex_results$species[to_idx])
    
    out <- furrr::future_pmap(params, function(from_idx, to_idx, E, from_var, to_var)
    {
        # pull out variables from the original block
        lib_ts <- simplex_results[[from_idx, "data"]]$abundance
        pred_ts <- simplex_results[[to_idx, "data"]]$abundance
        
        # compute CCM for actual connection
        ccm_actual <- ccm_func(cbind(lib_ts, pred_ts), E) %>%
            dplyr::mutate(data_type = "actual")
        # generate surrogates and compute CCM
        surr_ts <- simplex_results[[from_idx, "surrogate_data"]]
        
        ccm_surr <- purrr::map_dfr(seq(NCOL(surr_ts)),
                                   ~ccm_func(cbind(surr_ts[, .], pred_ts), E)) %>%
            dplyr::mutate(data_type = "surrogate")
        
        # combine outputs
        ccm_out <- dplyr::bind_rows(ccm_actual, ccm_surr) %>%
            dplyr::mutate(lib_column = from_var,
                          target_column = to_var)
        ccm_out$E <- E
        return(ccm_out)
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::select(lib_column, target_column, data_type, dplyr::everything()) %>%
        dplyr::mutate_at(c("lib_column", "target_column", "data_type"), as.factor)
    
}

#' @rdname compute_ccm
#' @description `build_ccm_plan` encapsulates the calculations of `compute_ccm`
#'   in a \code{\link[drake]{drake}} plan. This seems to run a bit faster than 
#'   the `compute_ccm` function.
#' @return `build_ccm_plan` returns a \code{\link[drake]{drake}} plan with 
#'   targets for the helper function, the params to be used for ccm, and the 
#'   ccm_results (as from \code{\link{compute_ccm}}).
#' 
#' @export
build_ccm_plan <- function(lib_sizes = seq(10, 100, by = 10),
                           random_libs = TRUE, num_samples = 100,
                           replace = TRUE, RNGseed = 42,
                           silent = TRUE)
{
    drake::drake_plan(
        ccm_func = function(df, E) {
            rEDM::ccm(df, E = E, lib_sizes = !!lib_sizes,
                      random_libs = !!random_libs, num_samples = !!num_samples,
                      replace = !!replace,
                      lib_column = 1, target_column = 2,
                      RNGseed = !!RNGseed, silent = !!silent) %>%
                rEDM::ccm_means(na.rm = TRUE) %>%
                dplyr::select(lib_size, num_pred, rho, mae, rmse)},
        ccm_params = expand.grid(from_idx = seq(NROW(simplex_results)),
                                 to_idx = seq(NROW(simplex_results))) %>%
            dplyr::mutate(E = simplex_results$best_E[from_idx],
                          from_var = simplex_results$species[from_idx],
                          to_var = simplex_results$species[to_idx]),
        ccm_results = furrr::future_pmap(ccm_params,
                                         function(from_idx, to_idx, E, from_var, to_var) {
                                             # pull out variables from the original block
                                             lib_ts <- simplex_results[[from_idx, "data"]]$abundance
                                             pred_ts <- simplex_results[[to_idx, "data"]]$abundance
                                             
                                             # compute CCM for actual connection
                                             ccm_actual <- ccm_func(cbind(lib_ts, pred_ts), E) %>%
                                                 dplyr::mutate(data_type = "actual")
                                             # generate surrogates and compute CCM
                                             surr_ts <- simplex_results[[from_idx, "surrogate_data"]]
                                             
                                             ccm_surr <- purrr::map_dfr(seq(NCOL(surr_ts)),
                                                                        ~ccm_func(cbind(surr_ts[, .], pred_ts), E)) %>%
                                                 dplyr::mutate(data_type = "surrogate")
                                             
                                             # combine outputs
                                             ccm_out <- dplyr::bind_rows(ccm_actual, ccm_surr) %>%
                                                 dplyr::mutate(lib_column = from_var,
                                                               target_column = to_var)
                                             ccm_out$E <- E
                                             return(ccm_out)
                                         }) %>%
            dplyr::bind_rows() %>%
            dplyr::select(lib_column, target_column, data_type, dplyr::everything()) %>%
            dplyr::mutate_at(c("lib_column", "target_column", "data_type"), as.factor)
    )
}

#' @title Identify the significant CCM links
#' @description Using the output of \code{\link{compute_ccm}}, determine which 
#'   links are significant. Significant links (`x`` "causes" `y``) are where:
#'   \itemize{
#'     \item ccm rho for the actual time series `x` and `y`` is greater than the 
#'       `null_quantile` level of the surrogate data (at the largest library size)
#'     \item the increase in rho for the actual time series `x`` and `y`` between the 
#'       smallest and largest library sizes is greater than `delta_rho_threshold`
#'   }
#' @param ccm_results the output from `compute_ccm` that we generate the links 
#'   for
#' @param null_quantile the quantile of the surrogate which we desire the actual 
#'   ccm rho to exceed (the default value of `0.975` corresponds to the upper 
#'   cutoff of a 95\% CI)
#' @param delta_rho_threshold the absolute increase which we desire the actual 
#'   ccm rho values to exceed, when comparing the smallest and largest library 
#'   sizes
#' @return A tibble with the filtered significant links, along with the values 
#'   of the test statistics that were computed (`delta_rho` and 
#'   `rho_minus_upper_q_null`)
#'
#' @export
compute_ccm_links <- function(ccm_results,
                              null_quantile = 0.975, 
                              delta_rho_threshold = 0.1)
{
    # process the compiled ccm results and summarize by each link:
    #   delta_rho is the rho for actual data at largest library size -
    #                    rho at smallest library size
    #   rho_minus_upper_q_null computes the upper quantile rho for the surrogate 
    #     time series, then takes the difference between that value and the rho 
    #     for the actual data, then selects only the value at the largest 
    #     library size
    ccm_out <- ccm_results %>%
        dplyr::group_by(lib_column, target_column, E) %>%
        tidyr::nest() %>%
        dplyr::mutate(delta_rho = purrr::map_dbl(data, ~ 
                                                     dplyr::filter(., data_type == "actual") %>%
                                                     dplyr::summarize(dplyr::last(rho, order_by = lib_size) -
                                                                          dplyr::first(rho, order_by = lib_size)) %>%
                                                     as.numeric), 
                      rho_minus_upper_q_null = purrr::map_dbl(data, ~
                                                                  dplyr::group_by(., lib_size, data_type) %>%
                                                                  dplyr::summarize(upper_q = quantile(rho, null_quantile, na.rm = TRUE)) %>%
                                                                  tidyr::spread(data_type, upper_q) %>%
                                                                  dplyr::ungroup() %>% 
                                                                  dplyr::summarize(dplyr::last(actual - surrogate, order_by = lib_size)) %>%
                                                                  as.numeric)
        ) %>% 
        dplyr::arrange(lib_column) %>%
        dplyr::select(-data)
    
    # Filter links
    #   - must either be significant (passing thresholds)
    #   - or be from self to self (keeps univariate E for later analysis)
    
    dplyr::filter(ccm_out, 
                  (delta_rho > delta_rho_threshold & rho_minus_upper_q_null > 0) |
                      lib_column == target_column)
}