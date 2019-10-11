context("Check data simulations")

test_that("3-species food chain model", {
  sample_times <- seq(0, 500, by = 5)
  columns <- c("time", "x", "y", "z")

  expect_error(dat <- simulate_3sp_food_chain(sample_times = sample_times), NA)
  dat <- as.data.frame(dat)

  expect_equal(dim(dat), c(length(sample_times), length(columns)))
  expect_true(all(columns %in% names(dat)))
  expect_equal(dat$time, sample_times)
  attributes(dat) <- attributes(dat)[sort(names(attributes(dat)))]
  expect_known_hash(round(dat, 4), "4177476790")
})

test_that("Flour beetle model", {
  sample_times <- seq(0, 500, by = 5)
  columns <- c("time", "L", "P", "A")

  expect_error(dat <- simulate_LPA_flour_beetles(sample_times = sample_times), NA)
  dat <- as.data.frame(dat)

  expect_equal(dim(dat), c(length(sample_times), length(columns)))
  expect_true(all(columns %in% names(dat)))
  expect_equal(dat$time, sample_times)
  attributes(dat) <- attributes(dat)[sort(names(attributes(dat)))]
  expect_known_hash(round(dat, 4), "59f1a23358")
})

test_that("Resource competition model", {
  sample_times <- seq(0, 50, by = 1)
  columns <- c("time", "R1", "R2", "R3", "N1", "N2", "N3", "N4", "N5")

  expect_error(dat <- simulate_resource_competition(sample_times = sample_times), NA)
  dat <- as.data.frame(dat)

  expect_equal(dim(dat), c(length(sample_times), length(columns)))
  expect_true(all(columns %in% names(dat)))
  expect_equal(dat$time, sample_times)
  attributes(dat) <- attributes(dat)[sort(names(attributes(dat)))]
  expect_known_hash(round(dat, 2), "ef6b6da73f")
})
