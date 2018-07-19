library(tidyverse)
library(lubridate)
library(rEDM)
library(parallel)

#### function to fit periodic spline through data, using yearday as input x ----

make_surrogate_annual_spline <- function(x, y, num_surr = 100)
{
    # filter out NA first
    n <- length(x)
    idx <- which(is.finite(x) & is.finite(y))
    x <- x[idx]
    y <- y[idx]
    
    xx <- c(x - 365, x, x + 365)
    yy <- c(y, y, y)
    seasonal_F <- smooth.spline(xx, yy)
    seasonal_cyc <- predict(seasonal_F, x)$y
    seasonal_resid <- y - seasonal_cyc
    
    vals <- matrix(unlist(
        lapply(seq(num_surr), function(i) {
            seasonal_cyc + sample(seasonal_resid, n)
        })
    ), ncol = num_surr)
    
    out <- matrix(NA, nrow = n, ncol = num_surr)
    out[idx, ] <- vals 
}

#' @title compute_simplex
#' @description Run simplex projection models for each time series in the 
#'   `block` and save the output, to determine the best embedding dimension for 
#'   each column.
#' @param block A data.frame containing time series for the community. Each 
#'   column is a time series of abundances.
#' @param E The embedding dimension or range of embedding dimensions to search 
#'   over. (default is `1:10`)
#' @param ... Remaining arguments are passed directly into `rEDM::simplex()`
#' 
#' @return A tibble with columns for the species name (taken from the original 
#'   column names), the abundance time series for each species, the output from 
#'   `rEDM::simplex()` and the best embedding dimension, as determined by the E 
#'   that minimizes MAE
#' 
#' @export
compute_simplex <- function(block, E = 1:10)
{
    block %>%
        gather(species, abundance) %>%
        group_by(species) %>%
        nest() %>% 
        mutate(simplex_out = 
                   map(data, ~ simplex(.$abundance, E = E, silent = TRUE))) %>%
        mutate(best_E = map_int(simplex_out, ~ filter(., mae == min(mae)) %>%
                                    pull(E)))
}


#' @title compute_CCM
#' @description Run pairwise CCM starting from the simplex_output. simplex projection models for each time series in the 
#'   `block` and save the output, to determine the best embedding dimension for 
#'   each column.
#' @param block A data.frame containing time series for the community. Each 
#'   column is a time series of abundances.
#' @param E The embedding dimension or range of embedding dimensions to search 
#'   over. (default is `1:10`)
#' @param ... Remaining arguments are passed directly into `rEDM::simplex()`
#' 
#' @return A tibble with columns for the species name (taken from the original 
#'   column names), the abundance time series for each species, the output from 
#'   `rEDM::simplex()` and the best embedding dimension, as determined by the E 
#'   that minimizes MAE
#' 
#' @export
compute_CCM <- function(simplex_results, 
                        lib_sizes = seq(10, 100, by = 10), 
                        silent = TRUE, RNGseed = 42, 
                        ..., 
                        num_cores = 1)
{
    ccm_func <- function(block, E)
    {
        cbind(block) %>%
            ccm(E = E, lib_sizes = lib_sizes, 
                lib_column = 1, target_column = 2, 
                silent = silent, RNGseed = RNGseed, 
                ...) %>% 
            ccm_means(na.rm = TRUE) %>%
            select(lib_size, num_pred, rho, mae, rmse)
    }
    
    params <- expand.grid(from_idx = seq(NROW(simplex_results)), 
                          to_idx = seq(NROW(simplex_results)))
    params$E <- simplex_results$best_E[params$from_idx]
    params$from_var <- simplex_results$species[params$from_idx]
    params$to_var <- simplex_results$species[params$to_idx]
    
    out <- mclapply(seq(NROW(params)), function(i) {
        # get the param variables we want    
        from_idx <- params$from_idx[i]
        to_idx <- params$to_idx[i]
        E <- params$E[i]
        
        # pull out variables from the original block
        lib_ts <- simplex_results[[from_idx, "data"]]$abundance
        pred_ts <- simplex_results[[to_idx, "data"]]$abundance
        
        # compute CCM for actual connection
        ccm_actual <- ccm_func(cbind(lib_ts, pred_ts), E) %>%
            mutate(data = "actual")
        # generate surrogates and compute CCM
        surr_ts <- simplex_results[[from_idx, "surrogate_data"]]
        
        ccm_surr <- map_dfr(seq(NCOL(surr_ts)), 
                            ~ccm_func(cbind(surr_ts[, .], pred_ts), E)) %>%
            mutate(data = "surrogate")
        
        # combine outputs
        ccm_out <- bind_rows(ccm_actual, ccm_surr) %>%
            mutate(lib_column = params$from_var[i], 
                   target_column = params$to_var[i])
        ccm_out$E <- E
        return(ccm_out)
    }, mc.cores = num_cores)
    
    do.call(bind_rows, out) %>%
        select(lib_column, target_column, data, everything()) %>%
        mutate_at(c("lib_column", "target_column", "data", "lib_size"), as.factor)
}

# Define causal links as those for which:
# - the difference between actual rho and top of 95% CI (quantile = 0.975) for 
#   null surrogates at the largest library size is > 0 (indicates significance)
# - the difference between actual rho at the largest library size and that at 
#   the smallest library size is > 0.1
identify_ccm_links <- function(ccm_results, 
                               delta_rho_threshold = 0.1, 
                               null_quantile = 0.975)
{
    # process the compiled ccm results and summarize by each link:
    #   delta_rho is the rho for actual data at largest library size - 
    #                    rho at smallest library size
    #   surr_CI contains the upper quantile rho for the surrogate time series
    #     and the difference between that value and the actual rho at that 
    #     library size
    #   rho_minus_upper_q_null takes the rho_minus_upper_q value at the 
    #     largest library size
    ccm_out <- ccm_results %>%
        group_by(lib_column, target_column, E) %>%
        nest() %>%
        mutate(delta_rho = map_dbl(data, function(df) { # compute delta_rho
            df %>%
                filter(data == "actual") %>%
                summarize(last(rho, order_by = lib_size) - 
                              first(rho, order_by = lib_size)) %>%
                as.numeric()
        }), 
        surr_CI = map(data, function(df) { # compute upper_q 
            left_join(df %>%
                          filter(data == "surrogate") %>% 
                          group_by(lib_size) %>%
                          summarize(upper_q = quantile(rho, null_quantile)), 
                      df %>% 
                          filter(data == "actual") %>%
                          select(lib_size, rho), 
                      by = "lib_size") %>%
                mutate(rho_minus_upper_q = rho - upper_q)
        }), 
        rho_minus_upper_q_null = map_dbl(surr_CI, function(df) {
            last(df$rho_minus_upper_q, order_by = df$lib_size)
        })) %>%
        arrange(lib_column) %>% 
        select(-data, -surr_CI)
    
    # Filter links
    #   - must either be significant (passing thresholds)
    #   - or be from self to self (keeps univariate E for later analysis)
    
    ccm_out %>%
        filter((delta_rho > delta_rho_threshold & rho_minus_upper_q_null > 0) | 
                   lib_column == target_column)
}

#' @title linear rescaling function
#' @description rescale the input vector to mean = 0, variance = 1
#' 
#' @param x A numeric vector.
#' @param na.rm A logical (default = TRUE) that affects how NAs are dealt with
#' 
#' @return A rescaled vector the same length as `x`
#' 
#' @export
rescale <- function(x, na.rm = TRUE) {(x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)}

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
#' 
#' @return A list with the matrix smap-coefficients for each predictor variable
#'   identified in CCM (these are the affected variables). The names in the list
#'   and the column names of the matrices use the variable names in the block.
#' 
#' @export
compute_smap_coeffs <- function(block, ccm_links)
{
    # rescale all the variables
    block <- block %>%
        mutate_all(rescale)
    
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
    smap_coeffs <- map(effect_variables, function(effect_var) {
        links <- ccm_links %>% filter(lib_column == !!effect_var)
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
                                make_block(block[, effect_var, drop = FALSE], 
                                           max_lag = num_effect_lags)[, -c(1, 2), drop = FALSE])
        }
        # determine best theta
        theta_list <- c(seq(0, 1, by = 0.1), seq(1.5, 10, by = 0.5)) #seq(0, 10, by = 0.1)
        theta_test <- block_lnlp(temp_block, method = "s-map", tp = 1,
                                 theta  = theta_list, silent = TRUE)
        best_theta <- theta_test$theta[which.min(theta_test$mae)]
        
        # re-run to get s-map coefficients
        smap_out <- block_lnlp(temp_block, method = "s-map", tp = 1,
                               theta  = best_theta, silent = TRUE, 
                               save_smap_coefficients = TRUE)
        
        smap_coeff <- smap_out$smap_coefficients[[1]]
        names(smap_coeff) <- c(names(temp_block), "const")
        return(smap_coeff)
    })
    names(smap_coeffs) <- effect_variables
    return(smap_coeffs)
}

#' @title Compute eigenvalues from the S-map coefficients
#' @description Using the S-map coefficients, assemble the appropriate 
#'   Jacobian matrices for each time point, and then compute the eigenvalues. 
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
#'   This function computes J at each time step and then the corresponding 
#'   eigenvalues if all the values are present (e.g. not NA).
#' @param smap_coeffs A list of the S-map coefficients matrices (as returned 
#'   from \code{\link{compute_smap_coeffs}}
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
        filter(lib_column %in% vars) %>%
        pull(E) %>%
        max()
    jacobian_dim <- max_E * length(vars)
    
    # check number of time points
    ts_length <- unique(map_dbl(smap_coeffs, NROW))
    stopifnot(length(ts_length) == 1)
    
    # for each predicted variable, 
    #   convert the time-indexed data.frame of coefficients into the 
    #   time-indexed matrix of coefficients
    smap_matrix_rows <- map(vars, function(var_name) {
        # get s-map coefficients without const column
        smap_df <- smap_coeffs[[var_name]] %>%
            select(-const)
        
        # output matrix (each row is 1 time step, each col is a var x lag)
        out <- matrix(0, nrow = ts_length, ncol = jacobian_dim)
        
        # find the correct column locations for the s-map coefficients
        regex_splits <- str_match(colnames(smap_df), "(.+)_(\\d+)|(.+)")
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
    
    smap_matrices <- map(seq(ts_length), function(t) {
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

compute_eigenvalues <- function(smap_matrices)
{
    eigenvalues <- map(smap_matrices, function(J) {
        if (any(is.na(J)))
            return(NA)
        return(eigen(J)$values)
    })
    return(eigenvalues)
}

