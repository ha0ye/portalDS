context("Check data simulations")

test_that("3-species food chain model", {
    sample_times <- seq(0, 500, by = 5)
    
    expect_error(dat <- simulate_3sp_food_chain(sample_times = sample_times), NA)
    dat <- as.data.frame(dat)
    
    expect_equal(dim(dat), c(101, 4))
    expect_true(all(c("time", "x", "y", "z") %in% names(dat)))
    expect_equal(dat$time, sample_times)
    attributes(dat) <- attributes(dat)[sort(names(attributes(dat)))]
    expect_known_hash(round(dat, 4), "4177476790")
})

test_that("Flour beetle model", {
    sample_times <- seq(0, 500, by = 5)
    
    expect_error(dat <- simulate_LPA_flour_beetles(sample_times = sample_times), NA)
    dat <- as.data.frame(dat)
    
    expect_equal(dim(dat), c(101, 4))
    expect_true(all(c("time", "L", "P", "A") %in% names(dat)))
    expect_equal(dat$time, sample_times)
    attributes(dat) <- attributes(dat)[sort(names(attributes(dat)))]
    expect_known_hash(round(dat, 4), "59f1a23358")
})

test_that("Resource competition model", {
    sample_times <- seq(0, 200, by = 2)
    
    expect_error(dat <- simulate_resource_competition(sample_times = sample_times), NA)
    dat <- as.data.frame(dat)
    
    expect_equal(dim(dat), c(101, 9))
    expect_true(all(c("time", "R1", "R2", "R3", "N1", "N2", "N3", "N4", "N5") %in% names(dat)))
    expect_equal(dat$time, sample_times)
    attributes(dat) <- attributes(dat)[sort(names(attributes(dat)))]
    expect_known_hash(round(dat, 3), "4e19e286cc")
})
