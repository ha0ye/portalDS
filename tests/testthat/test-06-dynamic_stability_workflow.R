context("Check full dynamic stability workflow")

test_that("Compute dynamic stability using block_3sp example data", {
  set.seed(616382247)

  # setup data
  data("block_3sp", package = "rEDM")
  var_names <- c("x", "y", "z")
  block <- setNames(
    block_3sp[, c("time", "x_t", "y_t", "z_t")],
    c("time", var_names)
  )

  # simplex results
  expect_error(simplex_results <- compute_simplex(block,
                                                  E_list = 3:5,
                                                  num_surr = 4,
                                                  surrogate_method = "random_shuffle", 
                                                  id_var = "time"
  ), NA)
  
  # round simplex outputs for hash
  simplex_results$simplex_out <- lapply(simplex_results$simplex_out, round, 4)
  expect_known_hash(simplex_results, "b00009fc2d")
  
  # ccm results
  expect_error(ccm_results <- compute_ccm(simplex_results,
    num_samples = 10
  ), NA)
  expect_equal(dim(ccm_results), c(450, 9))
  expect_setequal(as.character(ccm_results$lib_column), var_names)
  expect_setequal(as.character(ccm_results$target_column), var_names)
  expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
  expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
  expect_setequal(ccm_results$lib_size, seq(10, 100, 10))
  expect_known_hash(
    dplyr::mutate_at(ccm_results, c("rho", "mae", "rmse"), round, 4),
    "40a8eb1668"
  )

  # ccm links
  expect_error(ccm_links <- compute_ccm_links(ccm_results), NA)
  expect_known_hash(ccm_links, "70474772ed")

  # s-map coefficients
  expect_error(smap_coeffs <- compute_smap_coeffs(block, ccm_links, id_var = "time"), NA)
  expect_known_hash(lapply(smap_coeffs, round, 4), "c86a8cb782")

  expect_error(smap_matrices <- compute_smap_matrices(smap_coeffs, ccm_links), NA)
  expect_known_hash(lapply(smap_matrices, round, 4), "f837843845")

  expect_error(eigen_decomp <- compute_eigen_decomp(smap_matrices), NA)
  expect_error(eigenvalues <- eigen_decomp$values, NA)
  expect_known_hash(lapply(eigenvalues, round, 4), "0710a023e3")

  expect_error(eigenvectors <- eigen_decomp$vectors, NA)
  expect_known_hash(which(is.na(eigenvectors)), "bf89101447")
  expect_error(ev <- eigenvectors[!is.na(eigenvectors)], NA)
  expect_known_hash(
    lapply(ev, function(x) {
      s <- ifelse(Re(x[1, ]) == 0, 1, sign(Re(x[1, ])))
      s <- matrix(s, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
      round(s * x, 4)
    }),
    "c26e2df6f9"
  )
})

test_that("Verify compute_dynamic_stability function", {
  set.seed(616382247)

  # setup data
  data("block_3sp", package = "rEDM")
  var_names <- c("x", "y", "z")
  block <- setNames(
    block_3sp[, c("time", "x_t", "y_t", "z_t")],
    c("time", var_names)
  )

  # simplex results
  expect_error(results <- compute_dynamic_stability(block,
                                                    id_var = "time", 
                                                    E_list = 3:5,
                                                    surrogate_method = "random_shuffle",
                                                    num_surr = 4,
                                                    num_samples = 10), 
               NA)
  expect_is(results, "list")
  expect_true(all(c("block", "simplex_results", "ccm_results",
                    "ccm_links", "smap_matrices", "eigenvalues",
                    "eigenvectors") %in% 
                    names(results)))
  
  simplex_results <- results$simplex_results
  ccm_results <- results$ccm_results
  ccm_links <- results$ccm_links
  smap_coeffs <- results$smap_coeffs
  smap_matrices <- results$smap_matrices
  eigenvalues <- results$eigenvalues
  eigenvectors <- results$eigenvectors
  
  simplex_results$simplex_out <- lapply(simplex_results$simplex_out, round, 4)
  expect_known_hash(simplex_results, "b00009fc2d")
  expect_known_hash(
    dplyr::mutate_at(ccm_results, c("rho", "mae", "rmse"), round, 4),
    "40a8eb1668"
  )
  expect_known_hash(ccm_links, "70474772ed")
  expect_known_hash(lapply(smap_coeffs, round, 4), "c86a8cb782")
  expect_known_hash(lapply(smap_matrices, round, 4), "f837843845")
  expect_known_hash(lapply(eigenvalues, round, 4), "0710a023e3")
  expect_known_hash(which(is.na(eigenvectors)), "bf89101447")
  expect_error(ev <- eigenvectors[!is.na(eigenvectors)], NA)
  expect_known_hash(
    lapply(ev, function(x) {
      s <- ifelse(Re(x[1, ]) == 0, 1, sign(Re(x[1, ])))
      s <- matrix(s, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
      round(s * x, 4)
    }),
    "c26e2df6f9"
  )
})
