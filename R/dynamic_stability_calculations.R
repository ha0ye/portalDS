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
#'   \code{\link[rEDM]{make_surrogate_data}}
#' @param num_surr number of surrogates to compute
#' @param surr_params a list of named optional arguments to be passed into the 
#'   surrogate data function
#' @return A tibble with columns for the species name (taken from the original
#'   column names), the abundance time series for each species, the output from
#'   `rEDM::simplex()`, the best embedding dimension, as determined by the E
#'   that minimizes MAE, and surrogate time series
#'
#' @export
compute_simplex <- function(block, E_list = 1:10, 
                            surrogate_method = "annual_spline", 
                            num_surr = 100, 
                            surr_params = NULL)
{
    simplex_results <- block %>%
        dplyr::select(-censusdate) %>%
        tidyr::gather(species, abundance) %>%
        dplyr::group_by(species) %>%
        tidyr::nest() %>%
        dplyr::mutate(simplex_out =
                          purrr::map(data, ~ rEDM::simplex(.$abundance, E = E_list, silent = TRUE))) %>%
        dplyr::mutate(best_E = purrr::map_int(simplex_out, ~ dplyr::filter(., mae == min(mae)) %>%
                                                  dplyr::pull(E) %>%
                                                  as.integer))
    
    surrogate_method <- tolower(surrogate_method)
    if (surrogate_method == "twin")
    {
        simplex_results$surrogate_data <- 
            purrr::pmap(dplyr::select(simplex_results, data, best_E),
                        ~do.call(rEDM::make_surrogate_twin, 
                                 c(list(ts = ..1, dim = ..2, num_surr = num_surr), 
                                   surr_params)
                        )
            )
    } else if (surrogate_method == "annual_spline") {
        day_of_year <- lubridate::yday(block$censusdate)
        simplex_results$surrogate_data <- 
            purrr::pmap(dplyr::select(simplex_results, data),
                        ~do.call(make_surrogate_annual_spline, 
                                 c(list(ts = ..1, num_surr = num_surr), 
                                   surr_params)
                        )
            )
    } else {
        simplex_results$surrogate_data <- 
            purrr::pmap(dplyr::select(simplex_results, data),
                        ~do.call(rEDM::make_surrogate_data, 
                                 c(list(ts = ..1, num_surr = num_surr, 
                                        method = surrogate_method), 
                                   surr_params)
                        )
            )
    }
    return(simplex_results)
}

#' @title Run CCM on each pair of time series
#' @description Run pairwise CCM based on the simplex_output - using the best
#'   embedding dimension there, along with the surrogate time series

#' @param simplex_results the output of \code{\link{compute_simplex}}
#' @param num_cores the number of cores to use for computation
#' @inheritParams rEDM::ccm
#' @return A tibble with columns for the variables that we use in CCM, the data
#'   type (whether it's the "actual" time series or "surrogate"), the library
#'   size, and the results from CCM
#'
#' @export
compute_ccm <- function(simplex_results,
                        lib_sizes = seq(10, 100, by = 10),
                        random_libs = TRUE, num_samples = 100,
                        replace = TRUE, RNGseed = 42,
                        silent = TRUE, num_cores = 1)
{
    ccm_func <- function(block, E)
    {
        cbind(block) %>%
            rEDM::ccm(E = E, lib_sizes = lib_sizes,
                      random_libs = random_libs, num_samples = num_samples,
                      replace = TRUE,
                      lib_column = 1, target_column = 2,
                      RNGseed = RNGseed, silent = silent) %>%
            rEDM::ccm_means(na.rm = TRUE) %>%
            dplyr::select(lib_size, num_pred, rho, mae, rmse)
    }

    params <- expand.grid(from_idx = seq(NROW(simplex_results)),
                          to_idx = seq(NROW(simplex_results)))
    params$E <- simplex_results$best_E[params$from_idx]
    params$from_var <- simplex_results$species[params$from_idx]
    params$to_var <- simplex_results$species[params$to_idx]

    out <- parallel::mclapply(seq(NROW(params)), function(i) {
        # get the param variables we want
        from_idx <- params$from_idx[i]
        to_idx <- params$to_idx[i]
        E <- params$E[i]

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
            dplyr::mutate(lib_column = params$from_var[i],
                          target_column = params$to_var[i])
        ccm_out$E <- E
        return(ccm_out)
    }, mc.cores = num_cores)

    dplyr::bind_rows(out) %>%
        dplyr::select(lib_column, target_column, data_type, dplyr::everything()) %>%
        dplyr::mutate_at(c("lib_column", "target_column", "data_type"), as.factor)
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

#' @title Compute S-map coefficients for a given target variable
#' @description This function is meant to be called from within 
#'   \code{\link{compute_smap_coeffs}}, which also pre-generates the block so 
#'   that the first variable is to be predicted, and the remaining columns are 
#'   the causal variables and lags of the predicted variable. This function 
#'   searches over the values of theta for the best fit (by lowest MAE), and then 
#'   returns the data.frame with the s-map coefficients
#' @param block the input data with time delays already generated
#' @inheritParams rEDM::block_lnlp
#'
#' @return the data.frame with the s-map coefficients
#' @noRd
get_smap_coefficients <- function(block, 
                                  lib = c(1, NROW(block)),
                                  pred = c(1, NROW(block)),
                                  theta = c(seq(0, 1, by = 0.1), seq(1.5, 10, by = 0.5)))
{
    # determine best theta
    theta_test <- rEDM::block_lnlp(block, lib = lib, pred = pred,
                                   method = "s-map", tp = 1,
                                   theta  = theta, silent = TRUE)
    best_theta <- theta_test$theta[which.min(theta_test$mae)]
    
    # re-run to get s-map coefficients
    smap_out <- rEDM::block_lnlp(block, lib = lib, pred = pred,
                                 method = "s-map", tp = 1,
                                 theta  = best_theta, silent = TRUE,
                                 save_smap_coefficients = TRUE)
    return(smap_out$smap_coefficients[[1]])
}

#' @title Compute S-map coefficients for a community
#' @description Compute S-map models for each time series in the `block` and
#'   save out the coefficients. The coefficients represent the local linear
#'   model and can be used to infer properties of the system dynamics.
#' @details Suppose that there are causal links as follows:
#'   `x --> y`
#'   `x --> z`
#'   where `-->` indicates "cross-maps to". Then the interpretation from CCM is
#'   that `x` is affected by causes `y` and `z`. Thus, the predictive model for
#'   `x` should include `y` and `z` as predictors.
#'
#'   The S-map model is then setup as
#'   x_{t+1} = F(x_t, y_t, z_t, x_{t-1}, x_{t-2}, ...)
#'   where the number of predictors is equal to the best embedding dimension for
#'   `x`.
#' @param block A data.frame containing time series for the community. Each
#'   column is a time series of abundances.
#' @param ccm_links A data.frame containing the significant causal links. Each
#'   row is a causal link. The columns are:
#'   \describe{
#'   \item{`xmap_from`}{the column index of the predictor variable in CCM}
#'   \item{`xmap_to`}{the column index of the predicted variable in CCM}
#'   \item{`best_E`}{the best embedding dimension for CCM}
#'   }
#' @param rescale A logical, indicating whether to rescale each time series
#' @param rolling_forecast A logical, indicating whether to make individual
#'   rolling forecasts for the second half of the time series.
#'
#' @return A list with the matrix smap-coefficients for each predictor variable
#'   identified in CCM (these are the affected variables). The names in the list
#'   and the column names of the matrices use the variable names in the block.
#'
#' @export
compute_smap_coeffs <- function(block, ccm_links, rescale = TRUE,
                                rolling_forecast = FALSE)
{
    # rescale all the variables
    block <- block %>%
        dplyr::select(-censusdate)
    if (rescale)
    {
        block <- dplyr::mutate_all(block, norm_rescale)
    }

    if (is.factor(ccm_links$lib_column)) {
        ccm_links$lib_column <- as.character(ccm_links$lib_column)
    }
    if (is.factor(ccm_links$target_column)) {
        ccm_links$target_column <- as.character(ccm_links$target_column)
    }

    # check if column names are present or are valid indices
    if (is.numeric(ccm_links$lib_column))
    {
        stopifnot(min(ccm_links$lib_column) >= 1,
                  max(ccm_links$lib_column) <= NCOL(block))
    } else {
        stopifnot(ccm_links$lib_column %in% colnames(block))
    }
    if (is.numeric(ccm_links$target_column))
    {
        stopifnot(min(ccm_links$target_column) >= 1,
                  max(ccm_links$target_column) <= NCOL(block))
    } else {
        stopifnot(ccm_links$target_column %in% colnames(block))
    }

    effect_variables <- union(unique(ccm_links$lib_column),
                              unique(ccm_links$target_column))
    # Compute s-map coefficients for each variable
    #
    # Model setup:
    #   for each `xmap_from` variable,
    #   use the correct number of lags:
    #     all causal variables with no lag
    #     populate remainder with lags of predicted variable
    smap_coeffs <- purrr::map(effect_variables, function(effect_var) {
        links <- ccm_links %>% dplyr::filter(lib_column == !!effect_var)
        stopifnot(length(unique(links$E)) == 1) # check for unique best_E
        stopifnot(effect_var %in% links$target_column) # check for self-interaction

        E <- links$E[1]
        causal_var <- links$target_column
        causal_var <- c(effect_var, setdiff(causal_var, effect_var)) # reorder so effect_var is first

        # create temp_block
        #   how many total lags of effect_var do we need?
        num_effect_lags <- E - length(causal_var) + 1
        if (num_effect_lags < 1)
        {
            warning("Embedding dimension for ", effect_var, " was set at ", E,
                    " but ", length(causal_var), " predictors identified.\n",
                    " Using ", length(causal_var), " for E instead.")
            E <- length(causal_var)
            num_effect_lags <- 0
        }

        #   any causal vars and
        #   pad with lags of effect_var (dropping time and 0 lag columns)
        if (num_effect_lags == 0)
        {
            temp_block <- block[, causal_var, drop = FALSE]
        } else {
            temp_block <- cbind(block[, causal_var, drop = FALSE],
                                rEDM::make_block(block[, effect_var, drop = FALSE],
                                                 max_lag = num_effect_lags)[, -c(1, 2), drop = FALSE])
        }

        if (!rolling_forecast)
        {
            smap_coeff <- get_smap_coefficients(temp_block)
        } else {
            n <- NROW(block)
            lib <- c(1, floor(n/2))
            # initialize smap_coeff
            smap_coeff <- get_smap_coefficients(temp_block,
                                                lib = lib, pred = lib)
            to_fill <- data.frame(matrix(NA, nrow = n - NROW(smap_coeff), ncol = NCOL(smap_coeff)))
            names(to_fill) <- names(smap_coeff)
            smap_coeff <- rbind(smap_coeff, to_fill)

            # loop and generate new smap_coeff at each row
            for (lib_end in seq(floor(n/2) + 1, n))
            {
                lib <- c(1, lib_end)
                temp_smap_coeff <- get_smap_ceofficients(temp_block,
                                                         lib = lib, pred = lib)
                smap_coeff[lib_end - 1, ] <- temp_smap_coeff[lib_end - 1, ]
            }
        }
        names(smap_coeff) <- c(names(temp_block), "const")
        return(smap_coeff)
    })
    names(smap_coeffs) <- effect_variables
    return(smap_coeffs)
}

#' @title Generate the matrices of S-map coefficients
#' @description Using the S-map coefficients, assemble the appropriate
#'   Jacobian matrices for each time point
#' @details See \code{\link{compute_smap_coeffs}} for details on the input data.
#'   Let the variables in the system be x^{i} with i = 1..N.
#'   For the S-map model predicting x^{i}_{t+1}, let the coefficient
#'   corresponding to variable x^{j} at lag tau be c^{tau}_{ij}.
#'   Then the Jacobian is the block matrix:
#'   J = [[ C^0 C^1 ... C^K]
#'        [ I   0   ... 0  ]
#'        [ ... ... ... 0  ]
#'        [ ... ... I   0  ]]
#'   where K is the maximum lag, and C^{tau} is the matrix formed by the values
#'   c^{tau}_{ij}. (Note that many of these values will be 0.)
#'
#'   This function computes J at each time step.
#' @param smap_coeffs A list of the S-map coefficients for each predictor
#'   variable (as returned from \code{\link{compute_smap_coeffs}}
#' @param ccm_links A data.frame containing the significant causal links. Each
#'   row is a causal link. The columns are:
#'   \describe{
#'   \item{`xmap_from`}{the column index of the predictor variable in CCM}
#'   \item{`xmap_to`}{the column index of the predicted variable in CCM}
#'   \item{`best_E`}{the best embedding dimension for CCM}
#'   }
#'
#' @return A list with the matrix of smap-coefficients at each time point
#'
#' @export
compute_smap_matrices <- function(smap_coeffs, ccm_links)
{
    # identify only connected nodes
    ccm_links$lib_column <- as.factor(ccm_links$lib_column)
    ccm_links$target_column <- as.factor(ccm_links$target_column)

    vars <- union(levels(ccm_links$lib_column)[tabulate(ccm_links$lib_column) > 1],
                  levels(ccm_links$target_column)[tabulate(ccm_links$target_column) > 1])

    # get maximum E value (so we know size of the Jacobian matrix)
    max_E <- ccm_links %>%
        dplyr::filter(lib_column %in% vars) %>%
        dplyr::pull(E) %>%
        max()
    jacobian_dim <- max_E * length(vars)

    # check number of time points
    ts_length <- unique(purrr::map_dbl(smap_coeffs, NROW))
    stopifnot(length(ts_length) == 1)

    # for each predicted variable,
    #   convert the time-indexed data.frame of coefficients into the
    #   time-indexed matrix of coefficients
    smap_matrix_rows <- purrr::map(vars, function(var_name) {
        # get s-map coefficients without const column
        smap_df <- smap_coeffs[[var_name]] %>%
            dplyr::select(-const)

        # output matrix (each row is 1 time step, each col is a var x lag)
        out <- matrix(0, nrow = ts_length, ncol = jacobian_dim)

        # find the correct column locations for the s-map coefficients
        regex_splits <- stringr::str_match(colnames(smap_df), "(.+)_(\\d+)|(.+)")
        regex_splits[is.na(regex_splits[,2]), 2] <- ""  # convert NA into ""
        regex_splits[is.na(regex_splits[,3]), 3] <- "0" # convert NA into 0
        regex_splits[is.na(regex_splits[,4]), 4] <- ""  # convert NA into ""
        col_lookup <- data.frame(var = paste0(regex_splits[, 2], regex_splits[, 4]),
                                 lag = as.numeric(regex_splits[, 3]))
        col_lookup$col_idx <- match(col_lookup$var, vars) +
            col_lookup$lag * length(vars)

        # check for NAs
        if (any(is.na(col_lookup$col_idx)))
            warning("Computing column locations and found NAs for variable ", colnames(smap_df)[1])

        # populate s-map coefficients into output matrix
        smap_mat <- as.matrix(smap_df)
        out[, col_lookup$col_idx] <- smap_mat

        colnames(out) <- paste0(rep(vars, max_E), "_", rep(seq(max_E), each = length(vars)) - 1)
        return(out)
    })

    # for each time step
    #   create the matrix of smap_coeffs, check for NAs
    #   otherwise, fill the rest of the matrix accordingly, and compute the largest eigenvalue

    smap_matrices <- purrr::map(seq(ts_length), function(t) {
        # initialize J
        J <- matrix(0, nrow = jacobian_dim, ncol = jacobian_dim)

        # fill in rows for each predicted variable
        for (i in seq(smap_matrix_rows))
            J[i, ] <- (smap_matrix_rows[[i]])[t, ]

        # if there are NA values, then return NA
        if (any(is.na(J)))
            return(NA)

        # fill in identity matrix for lag relationships
        I_row_idx <- length(vars) + seq(length(vars) * (max_E - 1))
        I_col_idx <- seq(length(vars) * (max_E - 1))
        J[cbind(I_row_idx, I_col_idx)] <- 1

        colnames(J) <- colnames(smap_matrix_rows[[1]])
        rownames(J) <- colnames(J)

        return(J)
    })

    return(smap_matrices)
}

#' @title Compute the eigen-decomposition of the smap matrices

#' @param smap_matrices A list with the matrix of smap-coefficients at each
#'   time point \code{\link{compute_smap_matrices}}
#' @return A list with two elements:
#'   \describe{
#'     \item{`values`}{a list of the eigenvalues (a vector) for each time point}
#'     \item{`vectors`}{a list of the eigenvectors (a matrix, each column is an
#'       eigenvector) for each time point}
#'   }
#'
#' @export
compute_eigen_decomp <- function(smap_matrices)
{
    eigen_decomp <- purrr::map(smap_matrices, function(J) {
        if (any(is.na(J)))
            return(c("values" = NA, "vectors" = NA))
        out <- eigen(J)
        rownames(out$vectors) <- rownames(J)
        return(out)
    })
    return(purrr::transpose(eigen_decomp, .names = c("values", "vectors")))
}

