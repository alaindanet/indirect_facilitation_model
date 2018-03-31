context("Test that neighbor context function works")

parms <- c(
  z = 4,
  del = 0.1,
  b   = 0.8,
  cn  = 0.2,
  cp  = 0.2,
  cnp = 0.05,
  cpn = 0.2,
  g   = 0.08,
  n   = 1,
  m   = 0.2,
  gamma1 = 0.08)

rhoini <- 0.8

state_var <- c(
  N = rhoini * .5,
  P = rhoini * .5,
  NP = .1,
  NN = .1,
  PP = .1
  )

arg <- c(N = .4, NE = .1, E = .2, z = 4, del = .1, b = .8, c = .2, gamma1 = .05)

test_that("Check neighbor argument", {
  expect_error(check_nbs("NO"), "nbs is badly defined")
  expect_error(check_z("4"), "z is badly defined")
  expect_error(check_z(6), "z is badly defined")
  })

test_that("Neighbor context return expected errors",
  expect_error(Ncolonize(.4, .1, .2, 4, .1, .8, .2, .05, nbs = "NO"),
    "nbs is badly defined")
  )

test_that("NE context return the right proba", {
  expect_error(NE_context(nbs = "NE", NE = .1, E = .2, z = 4))
  expect_error(PE_context(nbs = "NE", NE = .1, E = .2, z = 4))
  })
