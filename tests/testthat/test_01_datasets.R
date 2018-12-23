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
    expect_error(portal_block <- readRDS(here::here("data/portal_block.RDS")), NA)
    expect_true(exists("portal_block"))
    expect_equal(dim(portal_block), c(440, 22))
    expect_true("censusdate" %in% names(portal_block))
    expect_equal(digest::digest(portal_block), 
                 "098292880fbb46a6365a20e9ed257936")
})

test_that("Portal Block 50 is correct", {
    expect_error(portal_block_50 <- readRDS(here::here("data/portal_block_50.RDS")), NA)
    expect_true(exists("portal_block_50"))
    expect_equal(dim(portal_block_50), c(440, 8))
    expect_true("censusdate" %in% names(portal_block_50))
    expect_equal(digest::digest(portal_block_50), 
                 "3d1dc7cc82df66e14cd8956b00022f1e")
})