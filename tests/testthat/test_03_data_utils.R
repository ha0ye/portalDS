context("Check data transformation utility functions")

test_that("norm_rescale works as expected", {
    set.seed(42)
    x <- runif(100)
    expect_error(x_rescaled <- norm_rescale(x), NA)
    expect_error(x_alt <- scale(x), NA)
    expect_equal(x_rescaled, as.numeric(x_alt))
})

test_that("make_surrogate_annual_spline works without errors", {
    set.seed(42)
    yday <- sample(seq(365), 50)
    ts <- rnorm(50)
    expect_error(out <- make_surrogate_annual_spline(yday, ts), NA)
    expect_is(out, "matrix")
    expect_equal(dim(out), c(50, 100))    
})

test_that("make_surrogate_annual_spline works performs correct calculations", {
    set.seed(43)
    yday <- sample(seq(365), 50)
    f <- function(yday)
    {
        5 * sin(yday * 2 * pi / 365) + 
            4 * cos(yday * 2 * pi / 200)
    }
    ts <- f(yday) + rnorm(50)
    expect_error(out <- make_surrogate_annual_spline(yday, ts), NA)
    expect_is(out, "matrix")
    expect_equal(dim(out), c(50, 100))
    err <- out - matrix(f(yday), ncol = 100, nrow = 50)
    expect_true(all(abs(err) < 3))
})
