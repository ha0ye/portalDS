#' @title Compute S-map coefficients for a given target variable
#' @description This function is meant to be called from within
#'   [compute_smap_coeffs()], which also pre-generates the block so that the 
#'   first variable is to be predicted, and the remaining columns are the causal 
#'   variables and lags of the predicted variable. This function searches over 
#'   the values of theta for the best fit (by lowest MAE), and then returns the 
#'   data.frame with the s-map coefficients
#' @param block the input data with time delays already generated
#' @inheritParams rEDM::block_lnlp
#'
#' @return the data.frame with the s-map coefficients
#' @export
get_smap_coefficients <- function(block,
                                  lib = c(1, NROW(block)),
                                  pred = c(1, NROW(block)),
                                  theta = c(seq(0, 1, by = 0.1), seq(1.5, 10, by = 0.5))) {
  # determine best theta
  theta_test <- rEDM::block_lnlp(block,
                                 lib = lib, pred = pred,
                                 method = "s-map", tp = 1,
                                 theta = theta, silent = TRUE
  )
  best_theta <- theta_test$theta[which.min(theta_test$mae)]
  
  # re-run to get s-map coefficients
  smap_out <- rEDM::block_lnlp(block,
                               lib = lib, pred = pred,
                               method = "s-map", tp = 1,
                               theta = best_theta, silent = TRUE,
                               save_smap_coefficients = TRUE
  )
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
#' @inheritParams compute_simplex
#' @return A list with the matrix smap-coefficients for each predictor variable
#'   identified in CCM (these are the affected variables). The names in the list
#'   and the column names of the matrices use the variable names in the block.
#'
#' @export
compute_smap_coeffs <- function(block, ccm_links, 
                                rescale = TRUE,
                                rolling_forecast = FALSE, 
                                id_var = NULL)
{
  # rescale all the variables
  abundances <- block
  if (!is.null(id_var))
  {
    abundances <- block %>%
      dplyr::select_at(dplyr::vars(-id_var))
  }
  if (rescale)
  {
    abundances <- dplyr::mutate_all(abundances, norm_rescale)
  }
  
  if (is.factor(ccm_links$lib_column)) {
    ccm_links$lib_column <- as.character(ccm_links$lib_column)
  }
  if (is.factor(ccm_links$target_column)) {
    ccm_links$target_column <- as.character(ccm_links$target_column)
  }
  
  # check if column names are present or are valid indices
  if (is.numeric(ccm_links$lib_column)) {
    stopifnot(
      min(ccm_links$lib_column) >= 1,
      max(ccm_links$lib_column) <= NCOL(abundances)
    )
  } else {
    stopifnot(ccm_links$lib_column %in% colnames(abundances))
  }
  if (is.numeric(ccm_links$target_column)) {
    stopifnot(
      min(ccm_links$target_column) >= 1,
      max(ccm_links$target_column) <= NCOL(abundances)
    )
  } else {
    stopifnot(ccm_links$target_column %in% colnames(abundances))
  }
  
  effect_variables <- union(
    unique(ccm_links$lib_column),
    unique(ccm_links$target_column)
  )
  # Compute s-map coefficients for each variable
  #
  # Model setup:
  #   for each `xmap_from` variable,
  #   use the correct number of lags:
  #     all causal variables with no lag
  #     populate remainder with lags of predicted variable
  smap_coeffs <- purrr::map(effect_variables, function(effect_var) {
    links <- ccm_links %>% dplyr::filter(.data$lib_column == !!effect_var)
    stopifnot(length(unique(links$E)) == 1) # check for unique best_E
    stopifnot(effect_var %in% links$target_column) # check for self-interaction
    
    E <- links$E[1]
    causal_var <- links$target_column
    causal_var <- c(effect_var, setdiff(causal_var, effect_var)) # reorder so effect_var is first
    
    # create temp_block
    #   how many total lags of effect_var do we need?
    num_effect_lags <- E - length(causal_var) + 1
    if (num_effect_lags < 1) {
      warning(
        "Embedding dimension for ", effect_var, " was set at ", E,
        " but ", length(causal_var), " predictors identified.\n",
        " Using ", length(causal_var), " for E instead."
      )
      E <- length(causal_var)
      num_effect_lags <- 0
    }
    
    #   any causal vars and
    #   pad with lags of effect_var (dropping time and 0 lag columns)
    if (num_effect_lags == 0) {
      temp_block <- abundances[, causal_var, drop = FALSE]
    } else {
      temp_block <- cbind(
        abundances[, causal_var, drop = FALSE],
        rEDM::make_block(abundances[, effect_var, drop = FALSE],
                         max_lag = num_effect_lags
        )[, -c(1, 2), drop = FALSE]
      )
    }
    
    if (!rolling_forecast) {
      smap_coeff <- get_smap_coefficients(temp_block)
    } else {
      n <- NROW(abundances)
      lib <- c(1, floor(n / 2))
      # initialize smap_coeff
      smap_coeff <- get_smap_coefficients(temp_block,
                                          lib = lib, pred = lib
      )
      to_fill <- data.frame(matrix(NA, nrow = n - NROW(smap_coeff), ncol = NCOL(smap_coeff)))
      names(to_fill) <- names(smap_coeff)
      smap_coeff <- rbind(smap_coeff, to_fill)
      
      # loop and generate new smap_coeff at each row
      for (lib_end in seq(floor(n / 2) + 1, n))
      {
        lib <- c(1, lib_end)
        temp_smap_coeff <- get_smap_coefficients(temp_block,
                                                 lib = lib, pred = lib
        )
        smap_coeff[lib_end - 1, ] <- temp_smap_coeff[lib_end - 1, ]
      }
    }
    names(smap_coeff) <- c(names(temp_block), "const")
    return(smap_coeff)
  })
  
  # add index labels for each matrix in list
  if (!is.null(id_var) && id_var %in% names(block))
  {
    labels <- block[[id_var]]
    smap_coeffs <- lapply(smap_coeffs, function(smap_coeff) {
      stopifnot(length(labels) == NROW(smap_coeff))
      rownames(smap_coeff) <- labels
      return(smap_coeff)
    })
  }
  
  names(smap_coeffs) <- effect_variables
  return(smap_coeffs)
}

#' @title Generate the matrices of S-map coefficients
#' @description Using the S-map coefficients, assemble the appropriate
#'   Jacobian matrices for each time point
#' @details See [compute_smap_coeffs()] for details on the input data.
#'   Let the variables in the system be x^{i} with i = 1..N.
#'   For the S-map model predicting x^{i}_{t+1}, let the coefficient
#'   corresponding to variable x^{j} at lag tau be c^{tau}_{ij}.
#'   Then the Jacobian is the block matrix, J = 
#'   \tabular{rrrrr}{
#'     C^0 \tab C^1 \tab ... \tab C^(d-1) \tab C^d\cr
#'       I \tab   0 \tab ... \tab       0 \tab   0\cr
#'       0 \tab   I \tab ... \tab       0 \tab   0\cr
#'     ... \tab ... \tab ... \tab     ... \tab ...\cr
#'       0 \tab   0 \tab ... \tab       I \tab   0
#'    }
#'   where d is the maximum lag, and C^{tau} is the matrix formed by the values
#'   c^{tau}_{ij}. (Note that many of these values will be 0.)
#'
#'   This function computes J at each time step.
#' @param smap_coeffs A list of the S-map coefficients for each predictor
#'   variable (as returned from [compute_smap_coeffs()])
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
  vars <- identify_connected_nodes(ccm_links)
  
  # restrict smap_coeffs to `vars`, and drop the constant coeff.
  smap_coeffs <- smap_coeffs[vars] %>%
    purrr::map(~ dplyr::select(., -.data$const))
  
  # generate the expanded state variables
  state_vars <- get_state_vars(smap_coeffs)
  
  # identify time steps to generate full matrices for
  valid_rows <- lapply(smap_coeffs, function(df) {
    range(which(apply(df, 1, function(x) {all(is.finite(x))})))
  })
  row_range <- seq(max(vapply(valid_rows, min, 0)), 
                   min(vapply(valid_rows, max, 0)))
  num_rows <- min(vapply(smap_coeffs, NROW, 1))
  
  # generate "blank" smap matrices
  smap_matrices <- generate_blank_smap_matrices(state_vars, row_range, num_rows)
  
  # for each target variable, copy all the coefficients into the correct matrix
  for (row_var in names(smap_coeffs))
  {
    # identify target row and column indices
    row_idx <- match(row_var, state_vars$name)
    coeff_matrix <- smap_coeffs[[row_var]]
    coeff_matrix <- as.matrix(coeff_matrix)
    col_var <- colnames(coeff_matrix)
    col_idx <- match(col_var, state_vars$name)
    matrix_idx <- cbind(row_idx, col_idx)
    
    # for each time step, copy values into correct matrix
    for (t in row_range)
    {
      smap_row <- coeff_matrix[t, ]
      smap_matrices[[t]][matrix_idx] <- smap_row
    }
  }
  
  # propagate list labels from `smap_coeffs`
  n <- length(smap_matrices)
  stopifnot(all(vapply(smap_coeffs, NROW, 0) == n))
  smap_coeff_names <- vapply(smap_coeffs, rownames, rep.int("", n))
  stopifnot(all(vapply(seq_len(NCOL(smap_coeff_names)), 
                       function(j) {identical(smap_coeff_names[,1], smap_coeff_names[,j])}, 
                       FALSE)
  ))
  names(smap_matrices) <- smap_coeff_names[, 1]
  
  return(smap_matrices)
}

# identify only connected nodes
# (variables that give or receive influence to other variables)
identify_connected_nodes <- function(ccm_links)
{
  ccm_links$lib_column <- as.factor(ccm_links$lib_column)
  ccm_links$target_column <- as.factor(ccm_links$target_column)
  
  union(
    levels(ccm_links$lib_column)[tabulate(ccm_links$lib_column) > 1],
    levels(ccm_links$target_column)[tabulate(ccm_links$target_column) > 1]
  )
}

# get the state variables used in the smap_coeffs
get_state_vars <- function(smap_coeffs)
{
  # organize the expanded state variables for the interaction matrix
  state_vars <- data.frame(name = unique(do.call(c, lapply(smap_coeffs, names))), 
                           stringsAsFactors = FALSE)
  out <- cbind(state_vars, 
               stringr::str_split(state_vars$name, "_", simplify = TRUE) %>%
                 as.data.frame(stringsAsFactors = FALSE) %>%
                 stats::setNames(c("base_var", "lag"))) %>%
    dplyr::mutate_at("lag", as.numeric) %>%
    tidyr::replace_na(list(lag = 0)) %>%
    dplyr::arrange(.data$lag, .data$base_var)
  stopifnot(setequal(out %>% dplyr::filter(lag == 0) %>% dplyr::pull(name), 
                     names(smap_coeffs)))
  return(out)
}

# generate the list of "blank" smap matrices
generate_blank_smap_matrices <- function(state_vars, row_range, num_rows)
{
  num_vars <- NROW(state_vars)
  
  # identify lagged cross-links
  lag_cross_links <- state_vars %>%
    dplyr::mutate(lag_name = paste(.data$base_var, .data$lag + 1, sep = "_"), 
                  row_idx = match(.data$lag_name, .data$name), 
                  col_idx = dplyr::row_number()) %>%
    dplyr::filter(!is.na(.data$row_idx))
  
  # generate the template matrix
  full_var_names <- paste(state_vars$base_var, state_vars$lag, sep = "_")
  M <- matrix(0, nrow = num_vars, ncol = num_vars, 
              dimnames = list(full_var_names, full_var_names))
  M[cbind(lag_cross_links$row_idx, lag_cross_links$col_idx)] <- 1
  
  # create list of smap matrices
  out <- as.list(rep(NA, num_rows))
  out[row_range] <- rep(list(M), length(row_range))
  return(out)
}
