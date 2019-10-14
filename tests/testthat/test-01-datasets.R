context("Check structure of included datasets")

test_that("Maizuru Block is correct", {
  expect_error(data("maizuru_block"), NA)
  expect_true(exists("maizuru_block"))
  expect_equal(dim(maizuru_block), c(285, 16))
  expect_true("censusdate" %in% names(maizuru_block))
  expect_known_hash(maizuru_block, "103576bf8c")
})

test_that("Portal Block is correct", {
  data_path <- system.file("extdata", "portal_block.RDS",
    package = "portalDS", mustWork = TRUE
  )
  expect_error(portal_block <- readRDS(data_path), NA)
  expect_true(exists("portal_block"))
  expect_equal(dim(portal_block), c(440, 22))
  expect_true("censusdate" %in% names(portal_block))
  attributes(portal_block) <- attributes(portal_block)[sort(names(attributes(portal_block)))]
  expect_known_hash(portal_block, "2efd6a078ea")
})

test_that("Portal Block 50 is correct", {
  data_path <- system.file("extdata", "portal_block_50.RDS",
    package = "portalDS", mustWork = TRUE
  )
  expect_error(portal_block_50 <- readRDS(data_path), NA)
  expect_true(exists("portal_block_50"))
  expect_equal(dim(portal_block_50), c(440, 8))
  expect_true("censusdate" %in% names(portal_block_50))
  attributes(portal_block_50) <- attributes(portal_block_50)[sort(names(attributes(portal_block_50)))]
  expect_known_hash(portal_block_50, "07a411a2f0")
})
