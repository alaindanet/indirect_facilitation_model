context("plot")

test_that("disconstinuity are correctly identified",{
  set.seed(12)
  series <- c(rnorm(10, 0, .01), rnorm(10, 1, .01))
  expect_equivalent(identify_discontinuity(series),
    c(rep(1, 10), rep(2, 10))
    )
})
