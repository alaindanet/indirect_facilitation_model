context("running_methods")

library(tibble)
set.seed(123)

run <- tibble(
  time = seq(0, 200),
  N = runif(n = 201),
  P = runif(n = 201)
  )
bad_run <- run %>%
  dplyr::mutate(
    N = c(NaN, runif(n = 200)),
    P = c(-400, runif(n = 200))
    )
test_that("is_normal_run works well", {
  expect_false(is_run_normal(bad_run))
  expect_true(is_run_normal(run))
    })


my_calc <- dplyr::slice(run, (n() - 100) :n()) %>%
  tidyr::gather(species, rho, -time) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(avg = mean(rho)) %>%
  tidyr::spread(species, avg)

test_that("avg_runs return a good summary of the simulation", {
  summary_run <- avg_runs(run, cut_row = 100)
  expect_equal(summary_run$N, my_calc$N, tolerance = .001)
  expect_equal(summary_run$P, my_calc$P, tolerance = .001)
  expect_true(summary_run$status)
  # Bad run:
  summary_run <- avg_runs(bad_run)
  expect_false(summary_run$status)
    })

test_that("avg_run returns NA if simulations has stoped early", {
  short_run <- tibble(
    time = seq(0, 2),
    N = runif(n = 3),
    P = runif(n = 3)
    )
  expect_equal(length(which(is.na(avg_runs(short_run)))), 3)
  expect_true(avg_runs(short_run)$status)

  short_bad_run <- tibble(
    time = seq(0, 2),
    N = rep(NaN, n = 3),
    P = rep(NaN, n = 3)
    )
  expect_equal(length(which(is.na(avg_runs(short_bad_run)))), 3)
  expect_false(avg_runs(short_bad_run)$status)
  })
test_that("the good states are returned", {
  expect_match(def_state(nurse = .3, protegee = .3, sim_status = TRUE),
    "coexistence")
  expect_match(def_state(nurse = .3, protegee = .3, sim_status = FALSE),
    "warning")
  expect_match(def_state(nurse = 3.84e-06, protegee = .3, sim_status = TRUE),
    "protégée")
  expect_match(def_state(nurse = .001, protegee = .001, sim_status = TRUE),
    "extinct")
  expect_match(def_state(nurse = .01, protegee = .001, sim_status = TRUE),
    "nurse")
  })


