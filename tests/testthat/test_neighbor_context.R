context("Test that neighbor context function works")
arg <- c(N = .4, NE = .1, E = .2, P = .4, PE = .1, z = 4, del = .1, b = .8,
  c = .2, gamma1 = .05)
test_that("Check neighbor argument", {
  expect_error(check_nbs("NO"), "nbs is badly defined")
  expect_null(check_nbs(NULL))
  expect_null(check_nbs("P"))
  expect_null(check_nbs("N"))
  expect_error(check_nbs("E"), "nbs is badly defined")
  expect_error(check_z("4"), "z is badly defined")
  expect_error(check_z(6), "z is badly defined")
  })
test_that("Neighbor context return expected errors",
  expect_error(Ncolonize(.4, .1, .2, 4, .1, .8, .2, .05, nbs = "NO"),
    "nbs is badly defined")
  )
test_that("NE context return the right proba", {
  expect_error(NE_context(nbs = "NE", NE = .1, E = .2, z = 4))
  expect_equal(NE_context(nbs = "N", NE = .1, E = .2, z = 4),
    with(as.list(arg), (1 / z + (z - 1) / z * NE / E)))
  expect_equal(NE_context(nbs = "P", NE = .1, E = .2, z = 4),
    with(as.list(arg), ( (z - 1) / z * NE / E)))
  expect_equal(NE_context(nbs = NULL, NE = .1, E = .2, z = 4),
    with(as.list(arg), NE/E))
  })
test_that("PE context return the right proba", {
  expect_error(PE_context(nbs = "NE", PE = .1, E = .2, z = 4))
  expect_equal(PE_context(nbs = "P", PE = .1, E = .2, z = 4),
    with(as.list(arg), (1 / z + (z - 1) / z * PE / E)))
  expect_equal(PE_context(nbs = "N", PE = .1, E = .2, z = 4),
    with(as.list(arg), ( (z - 1) / z * PE / E)))
  expect_equal(PE_context(nbs = NULL, PE = .1, E = .2, z = 4),
    with(as.list(arg), PE/E))
  })
