#' @title Dynamic stability analysis for community time series
#' @description [compute_dynamic_stability()] runs the full dynamic stability
#'   analysis. The analysis has multiple steps:
#'   \enumerate{
#'     \item{run simplex projection on each time series to identify the optimal
#'       embedding dimension}
#'     \item{generate surrogate time series, assumming that the data have just a
#'       seasonal pattern}
#'     \item{run ccm on each pairwise interaction, including the surrogate data}
#'     \item{identify the significant interactions by comparing the CCM for the
#'       real time series against the calculations for the surrogate data}
#'     \item{run S-map models for each time series, using the appropriate number of
#'       lags, and including the important interacting variables}
#'     \item{extract out the s-map coefficients from the models and assemble
#'       matrices for the system}
#'     \item{perform eigen-decomposition on the s-map coefficient matrices}
#'     \item{perform singular-value-decomposition on the s-map coefficient matrices}
#'     \item{compute volume contraction from the s-map coefficient matrices}
#'     \item{compute total variance from the s-map coefficient matrices}
#'   }
#' @param block a data.frame containing time series for the community. Each
#'   column is a time series of abundances, and a `censusdate` column is used
#'   as the column containing the time value.
#' @param results_file the location of the results to be stored on disk.
#' @param max_E largest E to examine using simplex projection; this sets the
#'   default range for `E_list`, but any setting for `E_list` will override the
#'   value for `max_E`
#' @inheritParams compute_simplex
#' @inheritParams compute_ccm
#' @inheritParams compute_smap_coeffs
#'
#' @return a list with named components for the individual output objects
#' XXX
#'
#' @export
compute_dynamic_stability <- function(block,
                                      results_file = NULL,
                                      id_var = "censudate", 
                                      max_E = 16, E_list = seq(max_E),
                                      surrogate_method = "annual_spline",
                                      num_surr = 200, surr_params = list(),
                                      lib_sizes = seq(10, 100, by = 10),
                                      random_libs = TRUE, num_samples = 100,
                                      replace = TRUE, RNGseed = 42,
                                      silent = TRUE, rescale = TRUE,
                                      rolling_forecast = FALSE)
{
  if (is.null(results_file) || !file.exists(results_file)) {
    results <- list()
  } else {
    results <- readRDS(results_file)
  }

  # check if data has been stored yet
  if (is.null(results$block)) {
    results$block <- block
  }

  # check for simplex results, compute if missing
  if (is.null(results$simplex_results)) {
    results$simplex_results <- compute_simplex(block,
      E_list = E_list,
      surrogate_method = surrogate_method,
      num_surr = num_surr,
      surr_params = surr_params
    )
  }

  # check for ccm results, compute if missing
  if (is.null(results$ccm_results)) {
    results$ccm_results <- compute_ccm(results$simplex_results,
      lib_sizes = lib_sizes,
      random_libs = random_libs,
      num_samples = num_samples,
      replace = replace,
      RNGseed = RNGseed,
      silent = silent
    )
  }

  # check for ccm links, compute if missing
  if (is.null(results$ccm_links)) {
    results$ccm_links <- compute_ccm_links(results$ccm_results)
  }

  # check for smap matrices, compute if missing
  if (is.null(results$smap_matrices)) {
    results$smap_coeffs <- compute_smap_coeffs(results$block, results$ccm_links,
      rescale = rescale,
      rolling_forecast = rolling_forecast
    )
    results$smap_matrices <- compute_smap_matrices(
      results$smap_coeffs,
      results$ccm_links
    )

    # add date labels for each matrix in list
    if (id_var %in% names(results$block))
    {
      stopifnot(length(results$smap_matrices) == NROW(results$block))
      names(results$smap_matrices) <- results$block[[id_var]]
    }
  }

  # check for eigenvalues
  if (is.null(results$eigenvalues) || is.null(results$eigenvectors)) {
    eigen_decomp <- compute_eigen_decomp(results$smap_matrices)
    results$eigenvalues <- eigen_decomp$values
    results$eigenvectors <- eigen_decomp$vectors
  }
  
  # check for svd
  if (is.null(results$singular_values) || is.null(results$singular_vectors_left)) {
    results$svd_decomp <- compute_svd_decomp(results$smap_matrices)
  }

  # check for volume contraction and total variance
  if (is.null(results$volume_contraction) || is.null(results$total_variance)) {
    results$volume_contraction <- compute_volume_contraction(results$smap_matrices)
    results$total_variance <- compute_total_variance(results$smap_matrices)
  }
  
  if (!is.null(results_file)) {
    saveRDS(results, file = results_file)
  }
  return(invisible(results))
}

#' @rdname compute_dynamic_stability
#' @title Create a drake plan for dynamic stability analysis
#' @description [build_dynamic_stability_plan()] creates a drake plan for the
#'   dynamic stability analysis.

#' @inheritParams compute_dynamic_stability
#' @export
build_dynamic_stability_plan <- function(max_E = 16, E_list = seq(max_E),
                                         surrogate_method = "annual_spline",
                                         num_surr = 200, surr_params = list(),
                                         lib_sizes = seq(10, 100, by = 10),
                                         random_libs = TRUE, num_samples = 100,
                                         replace = TRUE, RNGseed = 42,
                                         silent = TRUE, rescale = TRUE,
                                         rolling_forecast = FALSE) {
  drake::drake_plan(
    simplex_results = compute_simplex(
      block = block,
      E_list = !!E_list,
      surrogate_method = !!surrogate_method,
      num_surr = !!num_surr,
      surr_params = !!surr_params
    ),
    ccm_results = compute_ccm(
      simplex_results = simplex_results,
      lib_sizes = !!lib_sizes,
      random_libs = !!random_libs,
      num_samples = !!num_samples,
      replace = !!replace,
      RNGseed = !!RNGseed,
      silent = !!silent
    ),
    ccm_links = compute_ccm_links(ccm_results),
    smap_coeffs = compute_smap_coeffs(block, ccm_links,
      rescale = !!rescale,
      rolling_forecast = !!rolling_forecast
    ),
    smap_matrices = compute_smap_matrices(smap_coeffs, ccm_links),
    eigen_decomp = compute_eigen_decomp(smap_matrices),
    eigenvalues = eigen_decomp$values,
    eigenvectors = eigen_decomp$vectors, 
    svd_decomp = compute_svd_decomp(smap_matrices), 
    volume_contraction = compute_volume_contraction(smap_matrices), 
    total_variance = compute_total_variance(smap_matrices)
  )
}
