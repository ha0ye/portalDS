context("Check data obtaining functions work")

portal_data_path <- tempdir()

test_that("portalr downloads data ok", {
    expect_error(portalr::download_observations(path = portal_data_path, 
                                                version = "1.66.0"), NA)
})

test_that("make_portal_block works with default values", {
    expect_error(portal_block <- make_portal_block(path = portal_data_path), NA)
    expect_true(exists("portal_block"))
    expect_equal(dim(portal_block), c(440, 22))
    expect_true("censusdate" %in% names(portal_block))
    attributes(portal_block) <- attributes(portal_block)[sort(names(attributes(portal_block)))]
    expect_known_hash(portal_block, "2efd6a078e")
})

test_that("make_portal_block works with filtering by 50% present", {
    expect_error(portal_block_50 <- make_portal_block(path = portal_data_path, 
                                                      filter_q = 0.50), NA)
    expect_true(exists("portal_block_50"))
    expect_equal(dim(portal_block_50), c(440, 8))
    expect_true("censusdate" %in% names(portal_block_50))
    attributes(portal_block_50) <- attributes(portal_block_50)[sort(names(attributes(portal_block_50)))]
    expect_known_hash(portal_block_50, "07a411a2f0")
})