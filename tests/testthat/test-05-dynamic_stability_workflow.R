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
    expect_error(simplex_out <- compute_simplex(block, 
                                                E_list = E_list, 
                                                num_surr = 4, 
                                                surrogate_method = "random_shuffle"), NA)
    expect_known_hash(simplex_out, "ccf03f798a")
    
    # ccm results
    expect_error(ccm_results <- compute_ccm(simplex_out, 
                                            num_samples = 10), NA)
    expect_equal(dim(ccm_results), c(450, 9))
    expect_setequal(as.character(ccm_results$lib_column), var_names)
    expect_setequal(as.character(ccm_results$target_column), var_names)
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(ccm_results$lib_size, seq(10, 100, 10))
    expect_known_hash(ccm_results, "9f04e182a0")
    
    # ccm links
    expect_error(ccm_links <- identify_ccm_links(ccm_results), NA)
    expect_known_hash(ccm_links, "8450590e4c")
    
    # s-map coefficients
    
    # eigenvectors
    
    
    # eigenvalues
})