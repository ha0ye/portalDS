context("Check dynamic stability analysis functions work")

data_path <- system.file("extdata", "portal_block.RDS",
                         package = "portalDS", mustWork = TRUE)
portal_block <- readRDS(data_path)

data(maizuru_block)

test_that("compute_simplex works as expected", {
    E_list <- 3:5
    expect_error(simplex_out <- compute_simplex(portal_block, 
                                                E_list = E_list, 
                                                num_surr = 4), NA)
    # check columns
    expect_setequal(c("species", "data", "simplex_out", "best_E", "surrogate_data"), names(simplex_out))
    
    # check species values
    expect_setequal(simplex_out$species, setdiff(names(portal_block), "censusdate"))
    
    # check random data value
    idx <- sample(NROW(simplex_out), 1)
    simplex_out_data <- simplex_out[[idx, "data"]][[1]]
    original_portal_data <- portal_block[[simplex_out$species[idx]]]
    expect_equal(simplex_out_data, original_portal_data)
    
    # check simplex_out
    simplex_result <- simplex_out[[idx, "simplex_out"]]
    expect_equal(simplex_result$E, E_list)
    simplex_results <- simplex_out$simplex_out
    expect_known_hash(lapply(simplex_results, round, 4), "cbd599464a")
    
    # check best_E
    best_E <- simplex_out[["best_E"]]
    expect_true(all(best_E >= min(E_list)))
    expect_true(all(best_E <= max(E_list)))
    
    # check surrogate_data
    expect_error(surr_data <- simplex_out[[idx, "surrogate_data"]], NA)
    expect_equal(dim(surr_data), c(NROW(portal_block), 4))
})

test_that("compute_ccm run works as expected", {
    data("block_3sp", package = "rEDM")
    block <- block_3sp[, c("time", "x_t", "y_t", "z_t")]
    var_names <- c("x", "y", "z")
    names(block) <- c("censusdate", var_names)
    E_list <- 3:5
    expect_error(simplex_out <- compute_simplex(block, 
                                                E_list = E_list, 
                                                num_surr = 4, 
                                                surrogate_method = "random_shuffle"), NA)
    
    expect_error(ccm_results <- compute_ccm(simplex_out, 
                                            num_samples = 10), NA)
    expect_equal(dim(ccm_results), c(450, 9))
    expect_setequal(as.character(ccm_results$lib_column), var_names)
    expect_setequal(as.character(ccm_results$target_column), var_names)
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(as.character(ccm_results$data_type), c("actual", "surrogate"))
    expect_setequal(ccm_results$lib_size, seq(10, 100, 10))
    expect_known_hash(ccm_results, "9f04e182a0")
})

test_that("identify_ccm_links works as expected", {
    data_path <- system.file("extdata", "maizuru_ccm_results.RDS",
                             package = "portalDS", mustWork = TRUE)
    maizuru_ccm_results <- readRDS(data_path)
    
    expect_error(ccm_links <- identify_ccm_links(maizuru_ccm_results), NA)
    expect_equal(dim(ccm_links), c(32, 5))
    species_names <- setdiff(names(maizuru_block), "censusdate")
    expect_setequal(species_names, as.character(ccm_links$lib_column))
    expect_setequal(species_names,  as.character(ccm_links$target_column))
    expect_known_hash(ccm_links, "dd99bafb30")
})

## check remaining dynamic stability workflow functions