context("Plotting helper functions")

# vdiffr::manage_cases(filter = "plotting")

test_that("check matrix value data munging", {
  values_list <- list(a = c(3, 4), 
                      b = complex(real = c(1, 1), imaginary = c(-1, 1)))
  
  expect_error(v_df <- extract_matrix_values(values_list), NA)
  expect_equal(names(v_df), c("censusdate", "value", "rank"))
  expect_equal(v_df$censusdate, c("a", "a", "b", "b"))
  expect_equal(v_df$value, c(3, 4, 1-1i, 1+1i))
  expect_equal(v_df$rank, c(1, 2, 1, 2))
})
