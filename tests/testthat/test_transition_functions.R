context("Test that colonization function works")

rhoini <- 0.8
arg <- c(
  z = 4,
  del = 0.1,
  b   = 0.8,
  c  = 0.2,
  g   = 0.08,
  m   = 0.2,
  gamma1 = 0.08,
  tau = 20,
  N = rhoini * .5,
  P = rhoini * .5,
  NP = .1,
  NN = .1,
  PP = .1
  )
arg <- with(as.list(arg), c(arg,
    NE = N - NN - NP,
    PE = P - PP - NP,
    E = 1 - N - P,
    n = 1 - exp(- tau * gamma1)
    ))
arg

test_that("Ncolonization returns the right proba", {
  expect_equal(with(as.list(arg),
      Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "N")),
    with(as.list(arg),
      (del * N + (1 / z + ( (z - 1) / z) * NE / E) * (1 - del)) *
      (b - c * (1 - E) - gamma1))
    )
  expect_equal(with(as.list(arg),
      Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = "P")),
    with(as.list(arg),
      (del * N + ( ( (z - 1) / z) * NE / E) * (1 - del)) *
      (b - c * (1 - E) - gamma1))
    )
  expect_equal(with(as.list(arg),
      Ncolonize(N, NE, E, z, del, b, c, gamma1, nbs = NULL)),
    with(as.list(arg),
      (del * N + NE/E * (1 - del)) *
      (b - c * (1 - E) - gamma1))
    )
  })
test_that("Pcolonization returns the right proba", {
  expect_equal(with(as.list(arg),
      Pcolonize(P, N, NE, PE, E, z, del, b, c, g, n, nbs = "P")),
    with(as.list(arg),
      (del * P + (1 / z + ( (z - 1) / z) * PE / E) * (1 - del)) *
      (b - c * (1 - E) - g * (1 - ( ( (z - 1) / z) * NE / E) * n))
    )
    )
  expect_equal(with(as.list(arg),
      Pcolonize(P, N, NE, PE, E, z, del, b, c, g, n, nbs = "N")),
    with(as.list(arg),
      (del * P + ( ( (z - 1) / z) * PE / E) * (1 - del)) *
      (b - c * (1 - E) - g * (1 - (1 / z + ( (z - 1) / z) * NE / E) * n))
    )
    )
  expect_equal(with(as.list(arg),
      Pcolonize(P, N, NE, PE, E, z, del, b, c, g, n, nbs = NULL)),
    with(as.list(arg),
      (del * P + PE / E * (1 - del)) *
      (b - c * (1 - E) - g * (1 - (NE/E) * n))
    )
    )
  })
test_that("die returns the right proba", {
  expect_equal(with(as.list(arg), die(m)),
    with(as.list(arg), m)
    )
  })

