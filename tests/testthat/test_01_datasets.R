context("Check structure of included datasets")

test_that("Maizuru Block is correct", {
    expect_error(data("maizuru_block"), NA)
    expect_true(exists("maizuru_block"))
    expect_equal(dim(maizuru_block), c(285, 16))
    expect_true("censusdate" %in% names(maizuru_block))
    expect_equal(digest::digest(maizuru_block), 
                 "103576bf8c7bbafe82fd3a42296cf354")
})

test_that("Portal Block is correct", {
    data_path <- system.file("extdata", "portal_block.RDS",
                             package = "portalDS", mustWork = TRUE)
    expect_error(portal_block <- readRDS(data_path), NA)
    expect_true(exists("portal_block"))
    expect_equal(dim(portal_block), c(440, 22))
    expect_true("censusdate" %in% names(portal_block))
    attributes(portal_block) <- attributes(portal_block)[sort(names(attributes(portal_block)))]
    expect_equal(digest::digest(portal_block), 
                 "2efd6a078ea67c7373a0d4ab7e6331b1")
})

test_that("Portal Block 50 is correct", {
    data_path <- system.file("extdata", "portal_block_50.RDS",
                             package = "portalDS", mustWork = TRUE)
    expect_error(portal_block_50 <- readRDS(data_path), NA)
    expect_true(exists("portal_block_50"))
    expect_equal(dim(portal_block_50), c(440, 8))
    expect_true("censusdate" %in% names(portal_block_50))
    attributes(portal_block_50) <- attributes(portal_block_50)[sort(names(attributes(portal_block_50)))]
    expect_equal(digest::digest(portal_block_50), 
                 "07a411a2f03dd1eacfbad7366550596b")
})