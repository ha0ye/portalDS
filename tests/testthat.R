library(testthat)
library(portalDS)

if ("sample.kind" %in% names(formals(RNGkind))) {
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
RNGversion("3.5.2")

test_check("portalDS")
