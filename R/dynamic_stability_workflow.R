#' @title compute_dynamic_stability
#' @description Full dynamic stability analysis for community time series. The
#'   analysis has multiple components:
#'   (1) run simplex projection on each time series to identify the optimal
#'       embedding dimension
#'   (2) generate surrogate time series, assumming that the data have just a
#'       seasonal pattern
#'   (3) run ccm on each pairwise interaction, including the surrogate data
#'   (4) identify the significant interactions by comparing the CCM for the real
#'       time series against the calculations for the surrogate data
#'   (5) run S-map models for each time series, using the appropriate number of
#'       lags, and including the important interacting variables
#'   (6) extract out the s-map coefficients from the models and assemble
#'       matrices for the system
#'   (7) perform eigen-decomposition on the s-map coefficient matrices
#' @param block a data.frame containing time series for the community. Each
#'   column is a time series of abundances.
#' @param results_file the location of the results to be stored on disk.
#' @param max_E largest E to examine using simplex projection; this sets the
#'   default range for `E_list`, but any setting for `E_list` will override the
#'   value for `max_E`
#' @inheritParams compute_simplex
#' @inheritParams compute_ccm
#' @inheritParams compute_smap_coeffs
#'
#' @return none
#'
#' @export
compute_dynamic_stability <- function(block,
                                      results_file = "output/portal_ds_results.RDS",
                                      max_E = 16, E_list = seq(max_E), 
                                      surrogate_method = "annual_spline", num_surr = 200, surr_params = list(), 
                                      lib_sizes = c(6, 12, 24, 40, 80, 140, 220, 320, NROW(block)), 
                                      num_samples = 200, num_cores = 2,
                                      rescale = TRUE, rolling_forecast = FALSE)
{
    if (!file.exists(results_file))
    {
        results <- list()
    } else {
        results <- readRDS(results_file)
    }
    
    # check if data has been stored yet
    if (is.null(results$block))
    {
        results$block <- block
    }
    
    # check for simplex results, compute if missing
    if (is.null(results$simplex_results))
    {
        results$simplex_results <- compute_simplex(block,
                                                   E_list,
                                                   surrogate_method,
                                                   num_surr, 
                                                   surr_params)
    }
    
    # check for ccm results, compute if missing
    if (is.null(results$ccm_results))
    {
        results$ccm_results <- compute_ccm(results$simplex_results,
                                           lib_sizes = lib_sizes,
                                           num_samples = num_samples,
                                           num_cores = num_cores)
    }
    
    # check for ccm links, comput if missing
    if (is.null(results$ccm_links))
    {
        results$ccm_links <- compute_ccm_links(results$ccm_results)
    }
    
    # check for smap matrices, compute if missing
    if (is.null(results$smap_matrices))
    {
        smap_coeffs <- compute_smap_coeffs(results$block, results$ccm_links, 
                                           rescale, rolling_forecast)
        results$smap_matrices <- compute_smap_matrices(smap_coeffs,
                                                       results$ccm_links)
        
        # add date labels for each matrix in list
        stopifnot(length(results$smap_matrices) == NROW(results$block))
        names(results$smap_matrices) <- results$block$censusdate
    }
    
    # check for eigenvalues
    if (is.null(results$eigenvalues) || is.null(results$eigenvectors))
    {
        eigen_decomp <- compute_eigen_decomp(results$smap_matrices)
        results$eigenvalues <- eigen_decomp$values
        results$eigenvectors <- eigen_decomp$vectors
    }
    
    saveRDS(results, file = results_file)
}


#' @title Create a drake plan for dynamic stability analysis
#' @description Create a drake plan to perform dynamic stability analysis on 
#'   community time series data:
#'   (1) run simplex projection on each time series to identify the optimal
#'       embedding dimension
#'   (2) generate surrogate time series, assumming that the data have just a
#'       seasonal pattern
#'   (3) run ccm on each pairwise interaction, including the surrogate data
#'   (4) identify the significant interactions by comparing the CCM for the real
#'       time series against the calculations for the surrogate data
#'   (5) run S-map models for each time series, using the appropriate number of
#'       lags, and including the important interacting variables
#'   (6) extract out the s-map coefficients from the models and assemble
#'       matrices for the system
#'   (7) perform eigen-decomposition on the s-map coefficient matrices
#' @param block a data.frame containing time series for the community. Each
#'   column is a time series of abundances.
#' @param max_E largest E to examine using simplex projection; this sets the
#'   default range for `E_list`, but any setting for `E_list` will override the
#'   value for `max_E`
#' @inheritParams compute_simplex
#' @inheritParams compute_ccm
#' @inheritParams compute_smap_coeffs
#' @export
build_dynamic_stability_plan <- function(block,
                                         max_E = 16, E_list = seq(max_E), 
                                         surrogate_method = "annual_spline", num_surr = 200, surr_params = list(), 
                                         lib_sizes = c(6, 12, 24, 40, 80, 140, 220, 320, NROW(block)), 
                                         num_samples = 200)
{
    drake::drake_plan(
        simplex_results = compute_simplex(block = !!block,
                                          E_list = !!E_list,
                                          surrogate_method = !!surrogate_method,
                                          num_surr = !!num_surr, 
                                          surr_params = !!surr_params), 
        ccm_results = compute_ccm(simplex_results = simplex_results,
                                  lib_sizes = !!lib_sizes,
                                  num_samples = !!num_samples)
    )
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
