context("Check full dynamic stability workflow")

test_that("check dynamic stability using block_3sp example data", {
    set.seed(42)
    
    # setup data
    data("block_3sp", package = "rEDM")
    block <- block_3sp[, c("time", "x_t", "y_t", "z_t")]
    var_names <- c("x", "y", "z")
    names(block) <- c("censusdate", var_names)
    
    # simplex results
    E_list <- 3:5
    expect_error(simplex_results <- compute_simplex(block, 
                                                    E_list = E_list, 
                                                    num_surr = 4, 
                                                    surrogate_method = "random_shuffle"), NA)
    # round simplex outputs for hash
    simplex_results$simplex_out <- lapply(simplex_results$simplex_out, round, 4)
    expect_known_hash(simplex_results, "667a448396")
    
    # ccm results
    expect_error(ccm_results <- compute_ccm(simplex_results, 
                                            num_samples = 10), NA)
    expect_equal(dim(ccm_results), c(450, 9))
    expect_setequal(as.character(ccm_results$lib_column), var_names)
    expect_setequal(as.character(ccm_results$target_column), var_names)
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(ccm_results$lib_size, seq(10, 100, 10))
    expect_known_hash(dplyr::mutate_at(ccm_results, c("rho", "mae", "rmse"), round, 4), 
                      "4eb478eb8c")
    
    # ccm links
    expect_error(ccm_links <- compute_ccm_links(ccm_results), NA)
    expect_known_hash(ccm_links, "8450590e4c")
    
    # s-map coefficients
    expect_error(smap_coeffs <- compute_smap_coeffs(block, ccm_links), NA)
    expect_known_hash(lapply(smap_coeffs, round, 4), "3f88e8ddf9")
    
    expect_error(smap_matrices <- compute_smap_matrices(smap_coeffs, ccm_links), NA)
    smap_matrices <- lapply(smap_matrices, round, 4)
    expect_known_hash(smap_matrices, "5393993b55")
    
    expect_error(eigen_decomp <- compute_eigen_decomp(smap_matrices), NA)
    expect_error(eigenvalues <- eigen_decomp$values, NA)
    expect_known_hash(lapply(eigenvalues, round, 4), "05cbcc56e3")
    expect_error(eigenvectors <- eigen_decomp$vectors, NA)
    expect_known_hash(lapply(eigenvectors, round, 4), "43805666f4")
})