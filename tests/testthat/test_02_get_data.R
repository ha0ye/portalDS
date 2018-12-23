context("Check get_data produces identical results")

test_that("make_portal_block works with default values", {
    expect_error(portal_block <- make_portal_block(), NA)
    expect_true(exists("portal_block"))
    expect_equal(dim(portal_block), c(440, 22))
    expect_true("censusdate" %in% names(portal_block))
    expect_equal(digest::digest(portal_block), 
                 "098292880fbb46a6365a20e9ed257936")
})

test_that("make_portal_block works with filtering by 50% present", {
    expect_error(portal_block_50 <- make_portal_block(filter_q = 0.50), NA)
    expect_true(exists("portal_block_50"))
    expect_equal(dim(portal_block_50), c(440, 8))
    expect_true("censusdate" %in% names(portal_block_50))
    expect_equal(digest::digest(portal_block_50), 
                 "3d1dc7cc82df66e14cd8956b00022f1e")
})